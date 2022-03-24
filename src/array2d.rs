#[derive(Clone)]
pub struct Array3D<T> {
    width: usize,
    height: usize,
    length: usize,
    data: Vec<T>,
}

impl<T> Array3D<T> {
    pub fn from_array(width: usize, height: usize, length: usize, data: Vec<T>) -> Self {
        assert_eq!(data.len(), width * height * length);
        Self { width, height, length, data }
    }

    pub fn new(width: usize, height: usize, length: usize) -> Self
    where
        T: Default + Copy,
    {
        Self {
            width,
            height,
            length,
            data: vec![T::default(); width * height * length],
        }
    }

    pub fn data(&self) -> &[T] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    fn calc_index(&self, (x, y, z): (usize, usize, usize)) -> usize {
        debug_assert!(x < self.width);
        debug_assert!(y < self.height);
        debug_assert!(z < self.length);
        x + (y * self.width) + z * (self.width * self.height)
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn length(&self) -> usize {
        self.length
    }
}

impl<T> std::ops::Index<(usize, usize, usize)> for Array3D<T> {
    type Output = T;
    fn index(&self, pos: (usize, usize, usize)) -> &T {
        &self.data[self.calc_index(pos)]
    }
}

impl<T> std::ops::IndexMut<(usize, usize, usize)> for Array3D<T> {
    fn index_mut(&mut self, pos: (usize, usize, usize)) -> &mut T {
        let idx = self.calc_index(pos);
        &mut self.data[idx]
    }
}
