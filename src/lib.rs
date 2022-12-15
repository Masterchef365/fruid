pub type Array2D = idek_basics::Array2D<f32>;

pub struct FluidSim {
    u: Array2D,
    v: Array2D,
    smoke: Array2D,
}

impl FluidSim {
    pub fn new(width: usize, height: usize) -> Self {
        let arr = || Array2D::new(width, height);
        Self {
            u: Array2D::new(width + 1, height + 1),
            v: Array2D::new(width + 1, height + 1),
            smoke: Array2D::new(width, height),
        }
    }

    pub fn step(&mut self, dt: f32, overstep: f32) {
        todo!()
    }

    pub fn uv(&self) -> (&Array2D, &Array2D) {
        (&self.u, &self.v)
    }

    pub fn uv_mut(&mut self) -> (&mut Array2D, &mut Array2D) {
        (&mut self.u, &mut self.v)
    }

    pub fn smoke(&mut self) -> &mut Array2D {
        &mut self.smoke
    }

    pub fn width(&self) -> usize {
        self.u.width()
    }

    pub fn height(&self) -> usize {
        self.u.height()
    }
}
