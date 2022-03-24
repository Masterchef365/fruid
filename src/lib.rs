mod array2d;
pub type Array3D = array2d::Array3D<f32>;

// TODO: Transpose the traversal order; this one will be shit with regards to the cache

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
    /*
    let (nx, ny, nz) = inner_size(x);

    // Faces
    for i in 1..=ny {
        x[(0, i)] = if b == 1 { -x[(1, i)] } else { x[(1, i)] };
        x[(nx + 1, i)] = if b == 1 { -x[(nx, i)] } else { x[(nx, i)] };
    }

    // Faces
    for i in 1..=nx {
        x[(i, 0)] = if b == 2 { -x[(i, 1)] } else { x[(i, 1)] };
        x[(i, ny + 1)] = if b == 2 { -x[(i, ny)] } else { x[(i, ny)] };
    }

    // Corners
    x[(0, 0)] = 0.5 * (x[(1, 0)] + x[(0, 1)]);
    x[(0, ny + 1)] = 0.5 * (x[(1, ny + 1)] + x[(0, ny)]);
    x[(nx + 1, 0)] = 0.5 * (x[(nx, 0)] + x[(nx + 1, 1)]);
    x[(nx + 1, ny + 1)] = 0.5 * (x[(nx, ny + 1)] + x[(nx + 1, ny)]);
    */
    todo!()
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

fn lin_solve(b: i32, x: &mut Array3D, x0: &Array3D, scratch: &mut Array3D, a: f32, c: f32) {
    let (nx, ny, nz) = inner_size(x);

    let mut x = x;
    let mut scratch = scratch;

    for _ in 0..20 {
        for i in 1..=nx {
            for j in 1..=ny {
                for k in 1..=nz {
                    let neighbor_sum = neighbors((i, j, k)).iter().map(|&idx| x[idx]).sum::<f32>();
                    scratch[(i, j, k)] = (x0[(i, j, k)] + a * neighbor_sum) / c;
                }
            }
        }
        std::mem::swap(&mut scratch, &mut x);
        set_bnd(b, x);
    }
}

fn diffuse(b: i32, x: &mut Array3D, x0: &Array3D, scratch: &mut Array3D, diff: f32, dt: f32) {
    let (nx, ny, nz) = inner_size(x);
    let a = dt * diff * nx as f32 * ny as f32 * nz as f32;

    lin_solve(b, x, x0, scratch, a, 1. + 6. * a);
}

fn advect(b: i32, d: &mut Array3D, d0: &Array3D, u: &Array3D, v: &Array3D, w: &Array3D, dt: f32) {
    let (nx, ny, nz) = inner_size(d);

    let dt0 = dt * nx as f32;
    let dt1 = dt * ny as f32;
    let dt2 = dt * nz as f32;

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                let mut x = i as f32 - dt0 * u[(i, j, k)];
                let mut y = j as f32 - dt1 * v[(i, j, k)];
                let mut z = j as f32 - dt2 * w[(i, j, k)];

                if x < 0.5 {
                    x = 0.5
                };

                if x > nx as f32 + 0.5 {
                    x = nx as f32 + 0.5
                };

                let i0 = x as usize;
                let i1 = i0 + 1;

                if y < 0.5 {
                    y = 0.5
                };
                if y > ny as f32 + 0.5 {
                    y = ny as f32 + 0.5
                };

                let j0 = y as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f32;
                let s0 = 1. - s1;

                let t1 = y - j0 as f32;
                let t0 = 1. - t1;

                /*d[(i, j)] = 
                    + s0 * (t0 * d0[(i0, j0)] + t1 * d0[(i0, j1)])
                    + s1 * (t0 * d0[(i1, j0)] + t1 * d0[(i1, j1)]);
                */
                todo!()
            }
        }
    }
    set_bnd(b, d);
}

fn project(u: &mut Array3D, v: &mut Array3D, w: &mut Array3D, p: &mut Array3D, div: &mut Array3D, scratch: &mut Array3D) {
    let (nx, ny, nz) = inner_size(u);

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                let diffs = 
                    u[(i + 1, j, k)] - u[(i - 1, j, k)] 
                    + v[(i, j + 1, k)] - v[(i, j - 1, k)]
                    + w[(i, j, k + 1)] - w[(i, j, k - 1)];

                div[(i, j, k)] =
                    -0.5 * diffs / nx as f32; // TODO: Why is this just N? Shouldn't this be the total magnitude |nx, ny, nz|?

                p[(i, j, k)] = 0.0;
            }
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);

    lin_solve(0, p, div, scratch, 1., 6.);

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                u[(i, j, k)] -= 0.5 * nx as f32 * (p[(i + 1, j, k)] - p[(i - 1, j, k)]);
                v[(i, j, k)] -= 0.5 * ny as f32 * (p[(i, j + 1, k)] - p[(i, j - 1, k)]);
                w[(i, j, k)] -= 0.5 * ny as f32 * (p[(i, j, k + 1)] - p[(i, j, k - 1)]);
            }
        }
    }

    set_bnd(1, u);
    set_bnd(2, v);
    set_bnd(3, w);
}

fn dens_step(x: &mut Array3D, x0: &mut Array3D, u: &Array3D, v: &Array3D, w: &Array3D, scratch: &mut Array3D, diff: f32, dt: f32) {
    add_source(x, x0, dt);
    diffuse(0, x0, x, scratch, diff, dt);
    advect(0, x, x0, u, v, w, dt);
}

fn vel_step(
    u: &mut Array3D,
    v: &mut Array3D,
    w: &mut Array3D,
    u0: &mut Array3D,
    v0: &mut Array3D,
    w0: &mut Array3D,
    scratch: &mut Array3D,
    visc: f32,
    dt: f32,
) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);
    add_source(w, w0, dt);

    diffuse(1, u0, u, scratch, visc, dt);
    diffuse(2, v0, v, scratch, visc, dt);
    diffuse(3, w0, w, scratch, visc, dt);

    project(u0, v0, w0, u, v, scratch);

    //std::mem::swap(u0, u);
    //std::mem::swap(v0, v);

    advect(1, u, u0, u0, v0, w0, dt);
    advect(2, v, v0, u0, v0, w0, dt);
    advect(3, w, w0, u0, v0, w0, dt);

    project(u, v, w, u0, v0, scratch);
}

pub struct DensitySim {
    dens: Array3D,
    dens_prev: Array3D,
    scratch: Array3D,
}

impl DensitySim {
    pub fn new(width: usize, height: usize, length: usize) -> Self {
        let arr = || Array3D::new(width, height, length);
        Self {
            dens: arr(),
            dens_prev: arr(),
            scratch: arr(),
        }
    }

    pub fn step(&mut self, (u, v, w): (&Array3D, &Array3D, &Array3D), dt: f32, diff: f32) {
        dens_step(
            &mut self.dens,
            &mut self.dens_prev,
            u,
            v,
            w,
            &mut self.scratch,
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
    scratch: Array3D,
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
            scratch: arr(),
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
            &mut self.scratch,
            visc,
            dt,
        );
    }

    pub fn uvw(&self) -> (&Array3D, &Array3D) {
        (&self.u, &self.v)
    }

    pub fn uvw_mut(&mut self) -> (&mut Array3D, &mut Array3D) {
        (&mut self.u_prev, &mut self.v_prev)
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
