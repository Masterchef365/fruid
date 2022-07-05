pub type Array2D = idek_basics::Array2D<f32>;

// TODO: Transpose the traversal order; this one will be shit with regards to the cache

fn add_source(x: &mut Array2D, s: &Array2D, dt: f32) {
    x.data_mut()
        .iter_mut()
        .zip(s.data())
        .for_each(|(x, s)| *x += *s * dt);
}

fn inner_size(x: &Array2D) -> (usize, usize) {
    debug_assert!(x.width() >= 2);
    debug_assert!(x.height() >= 2);

    (x.width() - 2, x.height() - 2)
}

#[repr(i32)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Bounds {
    Positive = 0,
    NegX = 1,
    NegY = 2,
}

fn set_bnd(b: Bounds, x: &mut Array2D) {
    let (nx, ny) = inner_size(x);

    for i in 1..=ny {
        x[(0, i)] = if b == Bounds::NegX { -x[(1, i)] } else { x[(1, i)] };
        x[(nx + 1, i)] = if b == Bounds::NegX { -x[(nx, i)] } else { x[(nx, i)] };
    }

    for i in 1..=nx {
        x[(i, 0)] = if b == Bounds::NegY { -x[(i, 1)] } else { x[(i, 1)] };
        x[(i, ny + 1)] = if b == Bounds::NegY { -x[(i, ny)] } else { x[(i, ny)] };
    }

    x[(0, 0)] = 0.5 * (x[(1, 0)] + x[(0, 1)]);
    x[(0, ny + 1)] = 0.5 * (x[(1, ny + 1)] + x[(0, ny)]);
    x[(nx + 1, 0)] = 0.5 * (x[(nx, 0)] + x[(nx + 1, 1)]);
    x[(nx + 1, ny + 1)] = 0.5 * (x[(nx, ny + 1)] + x[(nx + 1, ny)]);
}

fn lin_solve(b: Bounds, x: &mut Array2D, x0: &Array2D, scratch: &mut Array2D, a: f32, c: f32) {
    let (nx, ny) = inner_size(x);

    let mut x = x;
    let mut scratch = scratch;

    for _ in 0..20 {
        for i in 1..=nx {
            for j in 1..=ny {

                let left = x[(i - 1, j)];
                let right = x[(i + 1, j)];
                let up = x[(i, j - 1)];
                let down = x[(i, j + 1)];

                let sum = left + right + up + down;

                let center_prev = x0[(i, j)];

                scratch[(i, j)] = (center_prev + a * sum) / c;
            }
        }
        std::mem::swap(&mut scratch, &mut x);
        set_bnd(b, x);
    }
}

fn diffuse(b: Bounds, x: &mut Array2D, x0: &Array2D, scratch: &mut Array2D, diff: f32, dt: f32) {
    let (nx, ny) = inner_size(x);
    let a = dt * diff * nx as f32 * ny as f32;

    lin_solve(b, x, x0, scratch, a, 1. + 4. * a);
}

fn mix(a: f32, b: f32, t: f32) -> f32 {
    (1. - t) * a + t * b
}

fn advect(b: Bounds, d: &mut Array2D, d0: &Array2D, u: &Array2D, v: &Array2D, dt: f32) {
    let (nx, ny) = inner_size(d);

    let dt0 = dt * nx as f32;
    let dt1 = dt * ny as f32;

    for i in 1..=nx {
        for j in 1..=ny {
            let mut x = i as f32 - dt0 * u[(i, j)];
            let mut y = j as f32 - dt1 * v[(i, j)];

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
            let t1 = y - j0 as f32;

            d[(i, j)] = mix(
                mix(d0[(i0, j0)], d0[(i0, j1)], t1),
                mix(d0[(i1, j0)], d0[(i1, j1)], t1),
                s1,
            );
        }
    }
    set_bnd(b, d);
}

enum DiffSide {
    X, Y
}

fn diff_acc(r: &mut Array2D, d: &Array2D, scale: f32, side: DiffSide) {
    let (nx, ny) = inner_size(d);
    for i in 1..=nx {
        for j in 1..=ny {
            let diff = match side {
                DiffSide::X => d[(i + 1, j)] - d[(i - 1, j)],
                DiffSide::Y => d[(i, j + 1)] - d[(i, j - 1)],
            };

            r[(i, j)] += scale * diff;
        }
    }
}

fn project(u: &mut Array2D, v: &mut Array2D, p: &mut Array2D, div: &mut Array2D, scratch: &mut Array2D) {
    let (nx, ny) = inner_size(u);

    for i in 1..=nx {
        for j in 1..=ny {
            div[(i, j)] = 0.;
            p[(i, j)] = 0.;
        }
    }

    diff_acc(div, u, -0.5 / nx as f32, DiffSide::X);
    diff_acc(div, v, -0.5 / ny as f32, DiffSide::Y);

    set_bnd(Bounds::Positive, div);
    set_bnd(Bounds::Positive, p);

    lin_solve(Bounds::Positive, p, div, scratch, 1., 4.);

    diff_acc(u, p, -0.5 * nx as f32, DiffSide::X);
    diff_acc(v, p, -0.5 * ny as f32, DiffSide::Y);

    set_bnd(Bounds::NegX, u);
    set_bnd(Bounds::NegY, v);
}

fn dens_step(x: &mut Array2D, x0: &mut Array2D, u: &Array2D, v: &Array2D, scratch: &mut Array2D, diff: f32, dt: f32) {
    //std::mem::swap(x0, x);

    diffuse(Bounds::Positive, x0, x, scratch, diff, dt);
    //std::mem::swap(x0, x);

    advect(Bounds::Positive, x, x0, u, v, dt);
}

fn vel_step(
    u: &mut Array2D,
    v: &mut Array2D,
    u0: &mut Array2D,
    v0: &mut Array2D,
    scratch: &mut Array2D,
    visc: f32,
    dt: f32,
) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);

    //std::mem::swap(u0, u);
    diffuse(Bounds::NegX, u0, u, scratch, visc, dt);

    //std::mem::swap(v0, v);
    diffuse(Bounds::NegY, v0, v, scratch, visc, dt);

    project(u0, v0, u, v, scratch);

    //std::mem::swap(u0, u);
    //std::mem::swap(v0, v);

    advect(Bounds::NegX, u, u0, u0, v0, dt);
    advect(Bounds::NegY, v, v0, u0, v0, dt);

    project(u, v, u0, v0, scratch);
}

pub struct DensitySim {
    dens: Array2D,
    dens_prev: Array2D,
    scratch: Array2D,
}

impl DensitySim {
    pub fn from_grid(dens: Array2D) -> Self {
        let arr = || Array2D::new(dens.width(), dens.height());
        Self {
            dens_prev: arr(),
            scratch: arr(),
            dens,
        }
    }

    pub fn step(&mut self, (u, v): (&Array2D, &Array2D), dt: f32, diff: f32) {
        dens_step(
            &mut self.dens,
            &mut self.dens_prev,
            u,
            v,
            &mut self.scratch,
            diff,
            dt,
        );
    }

    pub fn density(&self) -> &Array2D {
        &self.dens
    }

    pub fn density_mut(&mut self) -> &mut Array2D {
        &mut self.dens_prev
    }

    pub fn solve_pde(&mut self, init: bool) {
        let tmp = self.dens.clone();
        solve_pde(
            &mut self.dens, 
            &self.dens_prev, 
            &mut self.scratch, 
            0.5, 
            init
        );
        self.dens_prev = tmp;
    }
}

pub struct FluidSim {
    u: Array2D,
    v: Array2D,
    u_prev: Array2D,
    v_prev: Array2D,
    scratch: Array2D,
}

impl FluidSim {
    pub fn new(width: usize, height: usize) -> Self {
        let arr = || Array2D::new(width, height);
        Self {
            u: arr(),
            v: arr(),
            u_prev: arr(),
            v_prev: arr(),
            scratch: arr(),
        }
    }

    pub fn step(&mut self, dt: f32, visc: f32) {
        vel_step(
            &mut self.u,
            &mut self.v,
            &mut self.u_prev,
            &mut self.v_prev,
            &mut self.scratch,
            visc,
            dt,
        );
    }

    pub fn uv(&self) -> (&Array2D, &Array2D) {
        (&self.u, &self.v)
    }

    pub fn uv_mut(&mut self) -> (&mut Array2D, &mut Array2D) {
        (&mut self.u_prev, &mut self.v_prev)
    }

    pub fn width(&self) -> usize {
        self.u.width()
    }

    pub fn height(&self) -> usize {
        self.u.height()
    }
}


// https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html
pub fn solve_pde(x: &mut Array2D, x0: &Array2D, scratch: &mut Array2D, courant: f32, init: bool) {
    let (nx, ny) = inner_size(&x);

    for i in 1..=nx {
        for j in 1..=ny {
            let sum = x[(i - 1, j)]
                + x[(i + 1, j)]
                + x[(i, j - 1)]
                + x[(i, j + 1)]
                - 4. * x[(i, j)];

            scratch[(i, j)] = -x0[(i, j)] + 2. * x[(i, j)] + 0.5 * courant * sum;
        }
    }
    std::mem::swap(x, scratch);
}


