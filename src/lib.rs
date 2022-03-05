mod array2d;
pub type Array2D = array2d::Array2D<f32>;

// TODO: Transpose the traversal order; this one will be shit with regards to the cache

fn add_source(x: &mut Array2D, s: &Array2D, dt: f32) {
    x.data_mut()
        .iter_mut()
        .zip(s.data())
        .for_each(|(x, s)| *x = *s * dt);
}

fn inner_size(x: &Array2D) -> (usize, usize) {
    debug_assert!(x.width() >= 2);
    debug_assert!(x.height() >= 2);

    (x.width() - 2, x.height() - 2)
}

fn set_bnd(b: i32, x: &mut Array2D) {
    let (nx, ny) = inner_size(x);

    for i in 1..=ny {
        x[(0, i)] = if b == 1 { -x[(1, i)] } else { x[(1, i)] };
        x[(nx + 1, i)] = if b == 1 { -x[(nx, i)] } else { x[(nx, i)] };
    }

    for i in 1..=nx {
        x[(i, 0)] = if b == 2 { -x[(i, 1)] } else { x[(i, 1)] };
        x[(i, ny + 1)] = if b == 2 { -x[(i, ny)] } else { x[(i, ny)] };
    }

    x[(0, 0)] = 0.5 * (x[(1, 0)] + x[(0, 1)]);
    x[(0, ny + 1)] = 0.5 * (x[(1, ny + 1)] + x[(0, ny)]);
    x[(nx + 1, 0)] = 0.5 * (x[(nx, 0)] + x[(nx + 1, 1)]);
    x[(nx + 1, ny + 1)] = 0.5 * (x[(nx, ny + 1)] + x[(nx + 1, ny)]);
}

fn lin_solve(b: i32, x: &mut Array2D, x0: &Array2D, a: f32, c: f32) {
    let (nx, ny) = inner_size(x);

    for _ in 0..20 {
        for i in 1..=nx {
            for j in 1..=ny {
                x[(i, j)] = (x0[(i, j)]
                    + a * (x[(i - 1, j)] + x[(i + 1, j)] + x[(i, j - 1)] + x[(i, j + 1)]))
                    / c;
            }
        }
        set_bnd(b, x);
    }
}

fn diffuse(b: i32, x: &mut Array2D, x0: &Array2D, diff: f32, dt: f32) {
    let (nx, ny) = inner_size(x);
    let a = dt * diff * nx as f32 * ny as f32;

    lin_solve(b, x, x0, a, 1. + 4. * a);
}

fn advect(b: i32, d: &mut Array2D, d0: &Array2D, u: &Array2D, v: &Array2D, dt: f32) {
    let (nx, ny) = inner_size(d);

    let dt0 = dt * nx as f32; // TODO: why is this just N and not N*N?

    for i in 1..=nx {
        for j in 1..=ny {
            let mut x = i as f32 - dt0 * u[(i, j)];
            let mut y = j as f32 - dt0 * v[(i, j)];

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

            d[(i, j)] = s0 * (t0 * d0[(i0, j0)] + t1 * d0[(i0, j1)])
                + s1 * (t0 * d0[(i1, j0)] + t1 * d0[(i1, j1)]);
        }
    }
    set_bnd(b, d);
}

fn project(u: &mut Array2D, v: &mut Array2D, p: &mut Array2D, div: &mut Array2D) {
    let (nx, ny) = inner_size(u);

    for i in 1..=nx {
        for j in 1..=ny {
            div[(i, j)] =
                -0.5 * (u[(i + 1, j)] - u[(i - 1, j)] + v[(i, j + 1)] - v[(i, j - 1)]) / nx as f32; // TODO: Why is this just N?
            p[(i, j)] = 0.0;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);

    lin_solve(0, p, div, 1., 4.);

    for i in 1..=nx {
        for j in 1..=ny {
            u[(i, j)] -= 0.5 * nx as f32 * (p[(i + 1, j)] - p[(i - 1, j)]);
            v[(i, j)] -= 0.5 * ny as f32 * (p[(i, j + 1)] - p[(i, j - 1)]);
        }
    }

    set_bnd(1, u);
    set_bnd(2, v);
}

fn dens_step(x: &mut Array2D, x0: &mut Array2D, u: &Array2D, v: &Array2D, diff: f32, dt: f32) {
    add_source(x, x0, dt);
    std::mem::swap(x0, x);

    diffuse(0, x, x0, diff, dt);
    std::mem::swap(x0, x);

    advect(0, x, x0, u, v, dt);
}

fn vel_step(
    u: &mut Array2D,
    v: &mut Array2D,
    u0: &mut Array2D,
    v0: &mut Array2D,
    visc: f32,
    dt: f32,
) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);

    std::mem::swap(u0, u);
    diffuse(1, u, u0, visc, dt);

    std::mem::swap(v0, v);
    diffuse(2, v, v0, visc, dt);

    project(u, v, u0, v0);

    std::mem::swap(u0, u);
    std::mem::swap(v0, v);

    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);

    project(u, v, u0, v0);
}

pub struct Simulation {
    u: Array2D,
    v: Array2D,
    u_prev: Array2D,
    v_prev: Array2D,
    dens: Array2D,
    dens_prev: Array2D,
}

impl Simulation {
    pub fn new(width: usize, height: usize) -> Self {
        let arr = || Array2D::new(width, height);
        Self {
            u: arr(),
            v: arr(),
            u_prev: arr(),
            v_prev: arr(),
            dens: arr(),
            dens_prev: arr(),
        }
    }

    pub fn step(&mut self, dt: f32, visc: f32, diff: f32) {
        vel_step(
            &mut self.u,
            &mut self.v,
            &mut self.u_prev,
            &mut self.v_prev,
            visc,
            dt,
        );

        dens_step(
            &mut self.dens,
            &mut self.dens_prev,
            &mut self.u,
            &mut self.v,
            diff,
            dt,
        );
    }

    /*
    pub fn dimensions(&self) -> (usize, usize) {
        (self.u.width(), self.u.height())
    }
    */

    pub fn uv(&self) -> (&Array2D, &Array2D) {
        (&self.u, &self.v)
    }

    pub fn density(&self) -> &Array2D {
        &self.dens
    }

    pub fn uv_mut(&mut self) -> (&mut Array2D, &mut Array2D) {
        (&mut self.u_prev, &mut self.v_prev)
    }

    pub fn density_mut(&mut self) -> &mut Array2D {
        &mut self.dens_prev
    }
}
