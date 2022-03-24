mod array2d;
pub type Array3D = array2d::Array3D<f32>;

use crossbeam::thread::ScopedJoinHandle;

// TODO: Transpose the traversal order; this one will be shit with regards to the cache

const LIN_SOLVE_STEPS: usize = 20;

fn add_source(x: &mut Array3D, s: &Array3D, dt: f32) {
    x.data_mut()
        .iter_mut()
        .zip(s.data())
        .for_each(|(x, s)| *x += *s * dt);
}

fn inner_size(x: &Array3D) -> (usize, usize, usize) {
    debug_assert!(x.width() >= 2);
    debug_assert!(x.height() >= 2);
    debug_assert!(x.length() >= 2);

    (x.width() - 2, x.height() - 2, x.length() - 2)
}

fn set_bnd(b: i32, x: &mut Array3D) {
    let (nx, ny, nz) = inner_size(x);

    // This is a combinatorics problem.

    // First we find all faces (i, j) and set them equal to the pixel above them, negated if their
    // dimension is equal to b

    // Second, we find all edges (i) and set them equal to the average of the two adjacent
    // (recently filled) faces

    // Third, we find all of the corners () and set them equal to the average of the three
    // adjacent faces

    type Coord = (usize, usize, usize);
    fn fill_face(
        x: &mut Array3D,
        dim_i: usize,
        dim_j: usize,
        coord_fn: fn(Coord) -> Coord,
        negative: bool,
    ) {
        for i in 1..=dim_i {
            for j in 1..=dim_j {
                x[coord_fn((0, i, j))] = if negative {
                    -x[coord_fn((1, i, j))]
                } else {
                    x[coord_fn((1, i, j))]
                };
                x[coord_fn((dim_i + 1, i, j))] = if negative {
                    -x[coord_fn((dim_i, i, j))]
                } else {
                    x[coord_fn((dim_i, i, j))]
                };
            }
        }
    }

    fill_face(x, ny, nz, |(d, i, j)| (d, i, j), b == 1);
    fill_face(x, nx, nz, |(d, i, j)| (i, d, j), b == 2);
    fill_face(x, nx, ny, |(d, i, j)| (i, j, d), b == 3);

    // TODO: Corners and edges!!!
    /*
    x[(0, 0, 0)] = (x[(1, 0, 0)] + x[(0, 1, 0)] + x[(0, 0, 1)]) / 3.;
    x[(0, ny + 1)] = (x[(1, ny + 1)] + x[(0, ny)]) / 3.;
    x[(nx + 1, 0)] = 0.5 * (x[(nx, 0)] + x[(nx + 1, 1)]);
    x[(nx + 1, ny + 1)] = 0.5 * (x[(nx, ny + 1)] + x[(nx + 1, ny)]);
    */
}

fn neighbors((i, j, k): (usize, usize, usize)) -> [(usize, usize, usize); 6] {
    [
        (i - 1, j, k),
        (i + 1, j, k),
        (i, j - 1, k),
        (i, j + 1, k),
        (i, j, k - 1),
        (i, j, k + 1),
    ]
}

fn diffuse(b: i32, x: &mut Array3D, x0: &Array3D, accel: &mut LinSolveAccel, diff: f32, dt: f32) {
    let (nx, ny, nz) = inner_size(x);
    let a = dt * diff * nx as f32 * ny as f32 * nz as f32;

    accel.lin_solve(b, x, x0, a, 1. + 6. * a);
}

fn mix(a: f32, b: f32, t: f32) -> f32 {
    a * (1. - t) + b * t
}

fn advect(b: i32, d: &mut Array3D, d0: &Array3D, u: &Array3D, v: &Array3D, w: &Array3D, dt: f32) {
    let (nx, ny, nz) = inner_size(d);

    let dt0 = dt * nx as f32;
    let dt1 = dt * ny as f32;
    let dt2 = dt * nz as f32;

    fn advect_bounds(mut v: f32, dim: usize) -> (usize, usize, f32) {
        if v < 0.5 {
            v = 0.5
        }

        if v > dim as f32 + 0.5 {
            v = dim as f32 + 0.5
        }

        let i0 = v as usize;
        let i1 = i0 + 1;

        let s = v - i0 as f32;

        (i0, i1, s)
    }

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                let x = i as f32 - dt0 * u[(i, j, k)];
                let y = j as f32 - dt1 * v[(i, j, k)];
                let z = k as f32 - dt2 * w[(i, j, k)];

                let (i0, i1, s) = advect_bounds(x, nx);
                let (j0, j1, t) = advect_bounds(y, ny);
                let (k0, k1, g) = advect_bounds(z, nz);

                d[(i, j, k)] = mix(
                    mix(
                        mix(d0[(i0, j0, k0)], d0[(i1, j0, k0)], s),
                        mix(d0[(i0, j1, k0)], d0[(i1, j1, k0)], s),
                        t,
                    ),
                    mix(
                        mix(d0[(i0, j0, k1)], d0[(i1, j0, k1)], s),
                        mix(d0[(i0, j1, k1)], d0[(i1, j1, k1)], s),
                        t,
                    ),
                    g,
                );
            }
        }
    }
    set_bnd(b, d);
}

fn project(
    u: &mut Array3D,
    v: &mut Array3D,
    w: &mut Array3D,
    p: &mut Array3D,
    div: &mut Array3D,
    accel: &mut LinSolveAccel,
) {
    let (nx, ny, nz) = inner_size(u);

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                let diffs = u[(i + 1, j, k)] - u[(i - 1, j, k)] + v[(i, j + 1, k)]
                    - v[(i, j - 1, k)]
                    + w[(i, j, k + 1)]
                    - w[(i, j, k - 1)];

                div[(i, j, k)] = -0.5 * diffs / nx as f32; // TODO: Why is this just N? Shouldn't this be the total magnitude |nx, ny, nz|?

                p[(i, j, k)] = 0.0;
            }
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);

    accel.lin_solve(0, p, div, 1., 6.);

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                u[(i, j, k)] -= 0.5 * nx as f32 * (p[(i + 1, j, k)] - p[(i - 1, j, k)]);
                v[(i, j, k)] -= 0.5 * ny as f32 * (p[(i, j + 1, k)] - p[(i, j - 1, k)]);
                w[(i, j, k)] -= 0.5 * nz as f32 * (p[(i, j, k + 1)] - p[(i, j, k - 1)]);
            }
        }
    }

    set_bnd(1, u);
    set_bnd(2, v);
    set_bnd(3, w);
}

fn dens_step(
    x: &mut Array3D,
    x0: &mut Array3D,
    u: &Array3D,
    v: &Array3D,
    w: &Array3D,
    accel: &mut LinSolveAccel,
    diff: f32,
    dt: f32,
) {
    add_source(x, x0, dt);
    diffuse(0, x0, x, accel, diff, dt);
    advect(0, x, x0, u, v, w, dt);
}

fn vel_step(
    u: &mut Array3D,
    v: &mut Array3D,
    w: &mut Array3D,
    u0: &mut Array3D,
    v0: &mut Array3D,
    w0: &mut Array3D,
    accel: &mut LinSolveAccel,
    visc: f32,
    dt: f32,
) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);
    add_source(w, w0, dt);

    diffuse(1, u0, u, accel, visc, dt);
    diffuse(2, v0, v, accel, visc, dt);
    diffuse(3, w0, w, accel, visc, dt);

    project(u0, v0, w0, u, v, accel);

    //std::mem::swap(u0, u);
    //std::mem::swap(v0, v);

    advect(1, u, u0, u0, v0, w0, dt);
    advect(2, v, v0, u0, v0, w0, dt);
    advect(3, w, w0, u0, v0, w0, dt);

    project(u, v, w, u0, v0, accel);
}

pub struct DensitySim {
    dens: Array3D,
    dens_prev: Array3D,
    accel: LinSolveAccel,
}

impl DensitySim {
    pub fn new(width: usize, height: usize, length: usize) -> Self {
        let arr = || Array3D::new(width, height, length);
        Self {
            dens: arr(),
            dens_prev: arr(),
            accel: LinSolveAccel::new(width, height, length, LIN_SOLVE_STEPS),
        }
    }

    pub fn step(&mut self, (u, v, w): (&Array3D, &Array3D, &Array3D), dt: f32, diff: f32) {
        dens_step(
            &mut self.dens,
            &mut self.dens_prev,
            u,
            v,
            w,
            &mut self.accel,
            diff,
            dt,
        );
    }

    pub fn density(&self) -> &Array3D {
        &self.dens
    }

    pub fn density_mut(&mut self) -> &mut Array3D {
        &mut self.dens_prev
    }
}

pub struct FluidSim {
    u: Array3D,
    v: Array3D,
    w: Array3D,
    u_prev: Array3D,
    v_prev: Array3D,
    w_prev: Array3D,
    accel: LinSolveAccel,
}

impl FluidSim {
    pub fn new(width: usize, height: usize, length: usize) -> Self {
        let arr = || Array3D::new(width, height, length);
        Self {
            u: arr(),
            v: arr(),
            w: arr(),
            u_prev: arr(),
            v_prev: arr(),
            w_prev: arr(),
            accel: LinSolveAccel::new(width, height, length, LIN_SOLVE_STEPS),
        }
    }

    pub fn step(&mut self, dt: f32, visc: f32) {
        vel_step(
            &mut self.u,
            &mut self.v,
            &mut self.w,
            &mut self.u_prev,
            &mut self.v_prev,
            &mut self.w_prev,
            &mut self.accel,
            visc,
            dt,
        );
    }

    pub fn uvw(&self) -> (&Array3D, &Array3D, &Array3D) {
        (&self.u, &self.v, &self.w)
    }

    pub fn uvw_mut(&mut self) -> (&mut Array3D, &mut Array3D, &mut Array3D) {
        (&mut self.u_prev, &mut self.v_prev, &mut self.w_prev)
    }

    pub fn width(&self) -> usize {
        self.u.width()
    }

    pub fn height(&self) -> usize {
        self.u.height()
    }

    pub fn length(&self) -> usize {
        self.u.length()
    }
}

struct LinSolveAccel {
    /// A collection of scratch buffers, data buffers, and their assigned min and max z values
    bufs: Vec<AccelBuffer>,
}

struct AccelBuffer {
    /// Minimum Z in this area
    min_z: usize,
    /// Maximum  Z in this area
    max_z: usize,
    /// Minimum Z read from this area
    inner_min_z: usize,
    /// Maximum Z read from this area
    inner_max_z: usize,
    /// Buffer read from and written to
    x: Array3D,
    /// Scratch buffer for computations.
    /// Really we use a backbuffer/frontbuffer approach under the hood.
    scratch: Array3D,
    /// Number of steps to advance
    steps: usize,
}

impl LinSolveAccel {
    /// Creates a linear solver which will iterate the given number of steps of the jacobi method.
    pub fn new(width: usize, height: usize, length: usize, steps: usize) -> Self {
        let par = std::thread::available_parallelism()
            .map(|v| v.get())
            .unwrap_or(1);

        let area_length = length / par;
        let bufs: Vec<AccelBuffer> = (0..par)
            .map(|cpu| {
                let inner_min_z = cpu * area_length;
                let inner_max_z = inner_min_z + area_length;

                let min_z = inner_min_z.checked_sub(steps).unwrap_or(0);
                let max_z = (inner_max_z + steps).min(length);

                let x = Array3D::new(width, height, length);
                let scratch = Array3D::new(width, height, length);

                AccelBuffer {
                    min_z,
                    max_z,
                    inner_min_z,
                    inner_max_z,
                    x,
                    scratch,
                    steps,
                }
            })
            .collect();

        Self { bufs }
    }

    fn lin_solve(&mut self, b: i32, x: &mut Array3D, x0: &Array3D, a: f32, c: f32) {
        // TODO: Write buffers

        let x_ref: &Array3D = x;
        let bufs = &mut self.bufs;
        crossbeam::thread::scope(|s| {
            let threads: Vec<ScopedJoinHandle<AccelBuffer>> = bufs
                .drain(..)
                .map(|mut buf| {
                    s.spawn(move |_| {
                        for k in buf.min_z..buf.max_z {
                            for j in 0..x_ref.height() {
                                for i in 0..x_ref.width() {
                                    buf.x[(i, j, k)] = x_ref[(i, j, k)];
                                }
                            }
                        }

                        buf.solve(b, x0, a, c);
                        buf
                    })
                })
                .collect();

            *bufs = threads.into_iter().map(|t| t.join().unwrap()).collect();
        })
        .unwrap();
        drop(x_ref);
        drop(bufs);

        for buf in &mut self.bufs {
            for k in buf.inner_min_z..buf.inner_max_z {
                for j in 0..x.height() {
                    for i in 0..x.width() {
                        x[(i, j, k)] = buf.x[(i, j, k)];
                    }
                }
            }
        }

        // TODO: Read buffers
    }
}

impl AccelBuffer {
    pub fn solve(&mut self, b: i32, x0: &Array3D, a: f32, c: f32) {
        let (nx, ny, _) = inner_size(x0);

        for _ in 0..self.steps {
            for k in self.min_z + 1..self.max_z - 1 {
                for j in 1..=ny {
                    for i in 1..=nx {
                        let neighbor_sum = neighbors((i, j, k))
                            .iter()
                            .map(|&idx| self.x[idx])
                            .sum::<f32>();
                        self.scratch[(i, j, k)] = (x0[(i, j, k)] + a * neighbor_sum) / c;
                    }
                }
            }
            std::mem::swap(&mut self.scratch, &mut self.x);
            set_bnd(b, &mut self.x);
        }
    }
}
