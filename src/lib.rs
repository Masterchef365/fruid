use idek_basics::Array2D;
use num_complex::Complex32;

#[derive(Clone)]
pub struct FluidState {
    u: Array2D<f32>,
    v: Array2D<f32>,
    smoke: Array2D<Complex32>,
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
        // Copy values from last frame, in order to solve incompressibility
        self.write.u = self.read.u.clone();
        self.write.v = self.read.v.clone();

        // Force incompressibility
        for _ in 0..n_iters {
            for y in 1..self.write.v.height() - 2 {
                for x in 1..self.write.u.width() - 2 {
                    let dx = self.write.u[(x + 1, y)] - self.write.u[(x, y)];
                    let dy = self.write.v[(x, y + 1)] - self.write.v[(x, y)];

                    let d = overstep * (dx + dy) / 4.;

                    self.write.u[(x, y)] += d;
                    self.write.u[(x + 1, y)] -= d;

                    self.write.v[(x, y)] += d;
                    self.write.v[(x, y + 1)] -= d;
                }
            }
        }

        // Swap buffers such that the write buffer contains old data we intend to overwrite 
        // and the read buffer contains the fluid with incompressibility solved already
        std::mem::swap(&mut self.read.u, &mut self.write.u);
        std::mem::swap(&mut self.read.v, &mut self.write.v);

        // Advect velocity (u component)
        for y in 1..self.read.u.height() - 1 {
            for x in 1..self.read.u.width() - 1 {
                let u = self.read.u[(x, y)];
                let v = interp(&self.read.v, x as f32 - 0.5, y as f32 + 0.0);

                let px = x as f32 - u * dt;
                let py = y as f32 - v * dt;

                self.write.u[(x, y)] = interp(&self.read.u, px, py);
            }
        }

        // Advect velocity (v component)
        for y in 1..self.read.v.height() - 1 {
            for x in 1..self.read.v.width() - 1 {
                let u = interp(&self.read.u, x as f32 + 0.0, y as f32 - 0.5);
                let v = self.read.v[(x, y)];

                let px = x as f32 - u * dt;
                let py = y as f32 - v * dt;

                self.write.v[(x, y)] = interp(&self.read.v, px, py);
            }
        }

        // Swap the written buffers back into read again
        std::mem::swap(&mut self.read.u, &mut self.write.u);
        std::mem::swap(&mut self.read.v, &mut self.write.v);

        // Advect smoke
        for y in 1..self.read.v.height() - 2 {
            for x in 1..self.read.v.width() - 2 {
                let u = interp(&self.read.u, x as f32 + 0.5, y as f32 - 0.5);
                let v = interp(&self.read.v, x as f32 - 0.5, y as f32 + 0.5);

                let px = x as f32 - u * dt;
                let py = y as f32 - v * dt;

                self.write.smoke[(x, y)] = cinterp(&self.read.smoke, px, py);
            }
        }

        std::mem::swap(&mut self.read.smoke, &mut self.write.smoke);
    }

    pub fn uv(&self) -> (&Array2D<f32>, &Array2D<f32>) {
        (&self.read.u, &self.read.v)
    }

    pub fn uv_mut(&mut self) -> (&mut Array2D<f32>, &mut Array2D<f32>) {
        (&mut self.read.u, &mut self.read.v)
    }

    pub fn smoke_mut(&mut self) -> &mut Array2D<Complex32> {
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
fn clerp(a: Complex32, b: Complex32, t: f32) -> Complex32 {
    (1. - t) * a + t * b
}

/// Linear interpolation
fn lerp(a: f32, b: f32, t: f32) -> f32 {
    (1. - t) * a + t * b
}

/// Bilinear interpolation of the given grid at the given coordinates
fn interp(grid: &Array2D<f32>, x: f32, y: f32) -> f32 {
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


/// Bilinear interpolation of the given grid at the given coordinates
fn cinterp(grid: &Array2D<Complex32>, x: f32, y: f32) -> Complex32 {
    // Bounds enforcement. No panics!
    let tl_x = (x.floor() as isize).clamp(0, grid.width() as isize - 1) as usize;
    let tl_y = (y.floor() as isize).clamp(0, grid.height() as isize - 1) as usize;

    // Get corners
    let tl = grid[(tl_x, tl_y)];
    let tr = grid[(tl_x + 1, tl_y)];
    let bl = grid[(tl_x, tl_y + 1)];
    let br = grid[(tl_x + 1, tl_y + 1)];

    // Bilinear interpolation
    clerp(
        clerp(tl, tr, x.fract()), // Top row
        clerp(bl, br, x.fract()), // Bottom row
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
