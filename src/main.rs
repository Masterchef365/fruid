use fruid::{Array2D, Simulation};
use idek::{prelude::*, IndexBuffer};
mod graphics_builder;
use graphics_builder::GraphicsBuilder;

fn main() -> Result<()> {
    launch::<_, TriangleApp>(Settings::default().vr_if_any_args())
}

struct TriangleApp {
    line_verts: VertexBuffer,
    line_indices: IndexBuffer,
    line_gb: GraphicsBuilder,
    line_shader: Shader,

    tri_verts: VertexBuffer,
    tri_indices: IndexBuffer,
    tri_gb: GraphicsBuilder,

    sim: Simulation,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let sim = Simulation::new(150, 150);

        let mut line_gb = GraphicsBuilder::new();
        let mut tri_gb = GraphicsBuilder::new();

        //draw_density(&mut tri_gb, sim.density());
        draw_velocity_squares(&mut tri_gb, sim.uv());
        draw_density(&mut line_gb, sim.density());

        let line_verts = ctx.vertices(&line_gb.vertices, true)?;
        let line_indices = ctx.indices(&line_gb.indices, true)?;

        let tri_verts = ctx.vertices(&tri_gb.vertices, true)?;
        let tri_indices = ctx.indices(&tri_gb.indices, true)?;

        let line_shader = ctx.shader(
            DEFAULT_VERTEX_SHADER,
            DEFAULT_FRAGMENT_SHADER,
            Primitive::Lines,
        )?;

        Ok(Self {
            line_verts,
            line_indices,
            line_gb,
            line_shader,

            tri_verts,
            tri_indices,
            tri_gb,

            sim,
        })
    }

    fn frame(&mut self, ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        self.line_gb.clear();
        self.tri_gb.clear();

        //let (u, v) = self.sim.uv_mut();
        //draw_velocity_squares(&mut self.tri_gb, (u, v));

        //ctx.update_indices(&self.tri_gb.indices);

        let d = self.sim.density_mut();
        let center = (d.width() / 2, d.height() / 2);
        d[center] = 8.;

        let (u, v) = self.sim.uv_mut();
        let time = ctx.start_time().elapsed().as_secs_f32();
        let x = center.0 as f32 * (time.cos() + 1.);
        u[(x as usize, center.1)] = -10.;
        v[(x as usize, center.1)] = -10.;

        self.sim.step(0.1, 0., 0.);

        draw_density(&mut self.tri_gb, self.sim.density());
        //draw_velocity_squares(&mut self.tri_gb, self.sim.uv());
        ctx.update_vertices(self.tri_verts, &self.tri_gb.vertices)?;

        Ok(vec![
            /*DrawCmd::new(self.line_verts)
            .indices(self.line_indices)
            .shader(self.line_shader),*/
            DrawCmd::new(self.tri_verts).indices(self.tri_indices),
        ])
    }

    fn event(&mut self, ctx: &mut Context, platform: &mut Platform, event: Event) -> Result<()> {
        idek::simple_ortho_cam_ctx(ctx, platform);
        idek::close_when_asked(platform, &event);
        Ok(())
    }
}

fn draw_density(b: &mut GraphicsBuilder, density: &Array2D) {
    let cell_width = 2. / density.width() as f32;
    let cell_height = 2. / density.height() as f32;

    for i in 0..density.width() {
        let i_frac = (i as f32 / density.width() as f32) * 2. - 1.;
        for j in 0..density.height() {
            let j_frac = (j as f32 / density.height() as f32) * 2. - 1.;

            let density = density[(i, j)];
            let color = [density; 3];

            let mut push = |dx: f32, dy: f32| {
                let pos = [i_frac + dx, j_frac + dy, 0.];
                b.push_vertex(Vertex::new(pos, color))
            };

            let tl = push(0., 0.);
            let tr = push(cell_width, 0.);

            let bl = push(0., cell_height);
            let br = push(cell_width, cell_height);

            b.push_indices(&[bl, tr, tl, bl, br, tr]);
            //b.push_indices(&[tl, tr, bl, tr, br, bl]);
        }
    }
}


fn draw_velocity_squares(b: &mut GraphicsBuilder, (u, v): (&Array2D, &Array2D)) {
    let cell_width = 2. / u.width() as f32;
    let cell_height = 2. / u.height() as f32;

    for i in 0..u.width() {
        let i_frac = (i as f32 / u.width() as f32) * 2. - 1.;
        for j in 0..u.height() {
            let j_frac = (j as f32 / u.height() as f32) * 2. - 1.;

            let color = [u[(i, j)], v[(i, j)], 0.];

            let mut push = |dx: f32, dy: f32| {
                let pos = [i_frac + dx, j_frac + dy, 0.];
                b.push_vertex(Vertex::new(pos, color))
            };

            let tl = push(0., 0.);
            let tr = push(cell_width, 0.);

            let bl = push(0., cell_height);
            let br = push(cell_width, cell_height);

            b.push_indices(&[bl, tr, tl, bl, br, tr]);
            //b.push_indices(&[tl, tr, bl, tr, br, bl]);
        }
    }
}
