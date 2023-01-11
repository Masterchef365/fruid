pub type Array2D = idek_basics::Array2D<f32>;

#[derive(Clone)]
pub struct FluidState {
    u: Array2D,
    v: Array2D,
    smoke: Array2D,
}

pub struct FluidSim {
    read: FluidState,
    write: FluidState,
}

impl FluidSim {
    pub fn new(width: usize, height: usize) -> Self {
        let empty = FluidState {
            u: Array2D::new(width + 1, height),
            v: Array2D::new(width, height + 1),
            smoke: Array2D::new(width, height),
        };

        Self {
            read: empty.clone(),
            write: empty,
        }
    }

    pub fn step(&mut self, dt: f32, overstep: f32, n_iters: u32) {
        // Force incompressibility
        for _ in 0..n_iters {
            for y in 1..self.read.v.height() - 2 {
                for x in 1..self.read.u.width() - 2 {
                    let dx = self.read.u[(x + 1, y)] - self.read.u[(x, y)];
                    let dy = self.read.v[(x, y + 1)] - self.read.v[(x, y)];

                    let d = overstep * (dx + dy) / 4.;

                    self.read.u[(x, y)] += d;
                    self.read.u[(x + 1, y)] -= d;

                    self.read.v[(x, y)] += d;
                    self.read.v[(x, y + 1)] -= d;
                }
            }

            // Boundary conditions
            let h = self.read.v.height();
            for (y, smpl) in [(1, 2), (h - 2, h - 3)] {
                for i in 0..self.read.v.width() {
                    self.read.v[(i, y)] = -self.read.v[(i, smpl)];
                }
            }

            let w = self.read.u.width();
            for (x, smpl) in [(1, 2), (w - 2, w - 3)] {
                for i in 0..self.read.u.height() {
                    self.read.u[(x, i)] = -self.read.u[(smpl, i)];
                }
            }
        }

        fn advect(
            u: &Array2D,
            v: &Array2D,
            last_u: &Array2D,
            last_v: &Array2D,
            mut x: f32,
            mut y: f32,
            dt: f32,
        ) -> (f32, f32) {
            let n = 10;
            let dt = dt / n as f32;

            for i in 0..=n {
                let i = i as f32 / n as f32;
                let yy = y - 0.5;
                let xx = x - 0.5;
                let u = lerp(bilinear(&u, x, yy), bilinear(&last_u, x, yy), i);
                let v = lerp(bilinear(&v, xx, y), bilinear(&last_v, xx, y), i);
                x -= u * dt;
                y -= v * dt;
            }

            (x, y)
        }

        let old_u = self.write.u.clone();
        let old_v = self.write.v.clone();

        // Advect velocity (u component)
        for y in 1..self.read.u.height() - 1 {
            for x in 1..self.read.u.width() - 1 {
                let (px, py) = advect(&self.read.u, &self.read.v, &old_u, &old_v, x as f32, y as f32 + 0.5, dt);
                self.write.u[(x, y)] = bilinear(&self.read.u, px, py - 0.5);
            }
        }

        // Advect velocity (v component)
        for y in 1..self.read.v.height() - 1 {
            for x in 1..self.read.v.width() - 1 {
                let (px, py) = advect(&self.read.u, &self.read.v, &old_u, &old_v, x as f32 + 0.5, y as f32, dt);
                self.write.v[(x, y)] = bilinear(&self.read.v, px - 0.5, py);
            }
        }

        // Swap the written buffers back into read again
        std::mem::swap(&mut self.read.u, &mut self.write.u);
        std::mem::swap(&mut self.read.v, &mut self.write.v);

        // Advect smoke
        let mut mass = 0.;
        for y in 1..self.read.v.height() - 2 {
            for x in 1..self.read.v.width() - 2 {
                let (px, py) = advect(
                    &self.read.u,
                    &self.read.v,
                    &old_u,
                    &old_v,
                    x as f32 + 0.5,
                    y as f32 + 0.5,
                    dt,
                );
                self.write.smoke[(x, y)] = bilinear(&self.read.smoke, px - 0.5, py - 0.5);
                mass += self.write.smoke[(x, y)];
            }
        }
        dbg!(mass);

        std::mem::swap(&mut self.read.smoke, &mut self.write.smoke);
    }

    pub fn uv(&self) -> (&Array2D, &Array2D) {
        (&self.read.u, &self.read.v)
    }

    pub fn uv_mut(&mut self) -> (&mut Array2D, &mut Array2D) {
        (&mut self.read.u, &mut self.read.v)
    }

    pub fn smoke_mut(&mut self) -> &mut Array2D {
        &mut self.read.smoke
    }

    pub fn width(&self) -> usize {
        self.read.u.width()
    }

    pub fn height(&self) -> usize {
        self.read.u.height()
    }
}

/// Linear interpolation
fn lerp(a: f32, b: f32, t: f32) -> f32 {
    (1. - t) * a + t * b
}

/// Bilinear interpolation of the given grid at the given coordinates
#[track_caller]
fn bilinear(grid: &Array2D, x: f32, y: f32) -> f32 {
    // Bounds enforcement. No panics!
    let tl_x = (x.floor() as isize).clamp(0, grid.width() as isize - 1) as usize;
    let tl_y = (y.floor() as isize).clamp(0, grid.height() as isize - 1) as usize;

    // Get corners
    let tl = grid[(tl_x, tl_y)];
    let tr = grid[(tl_x + 1, tl_y)];
    let bl = grid[(tl_x, tl_y + 1)];
    let br = grid[(tl_x + 1, tl_y + 1)];

    // Bilinear interpolation
    lerp(
        lerp(tl, tr, x.fract()), // Top row
        lerp(bl, br, x.fract()), // Bottom row
        y.fract(),
    )
}

/*
   fn enforce_bounds() {
// Set grid boundaries
for x in 0..self.write.u.width() {
let top = (x, 0);
let bottom = (x, self.write.u.height() - 1);
self.write.u[top] = 0.;
self.write.u[bottom] = 0.;
self.write.v[top] = 0.;
self.write.v[bottom] = 0.;
}

for y in 0..self.write.u.height() {
let left = (0, y);
let right = (self.write.u.width() - 1, y);
self.write.u[left] = 0.;
self.write.u[right] = 0.;
self.write.v[left] = 0.;
self.write.v[right] = 0.;
}
}
*/
