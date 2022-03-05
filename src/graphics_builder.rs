use idek::Vertex;

#[derive(Default, Clone, Debug)]
pub struct GraphicsBuilder {
    pub vertices: Vec<Vertex>,
    pub indices: Vec<u32>,
}

impl GraphicsBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Push a Vertex and return it's index
    pub fn push_vertex(&mut self, v: Vertex) -> u32 {
        let idx: u32 = self
            .vertices
            .len()
            .try_into()
            .expect("Vertex limit exceeded");
        self.vertices.push(v);
        idx
    }

    /// Push an index
    pub fn push_indices(&mut self, idx: &[u32]) {
        self.indices.extend_from_slice(idx);
    }

    /// Erase all content
    pub fn clear(&mut self) {
        self.indices.clear();
        self.vertices.clear();
    }
}
