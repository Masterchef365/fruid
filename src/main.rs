use fruid::array2d::Array2D;
use idek::{prelude::*, IndexBuffer};
mod graphics_builder;
use graphics_builder::GraphicsBuilder;
use rand::prelude::*;

fn main() -> Result<()> {
    launch::<_, TriangleApp>(Settings::default().vr_if_any_args())
}

const DENSITY_Z: f32 = 0.5;
const VELOCITY_Z: f32 = 0.0;

struct TriangleApp {
    //line_verts: VertexBuffer,
    //line_indices: IndexBuffer,
    //line_gb: GraphicsBuilder,
    //line_shader: Shader,

    tri_verts: VertexBuffer,
    tri_indices: IndexBuffer,
    tri_gb: GraphicsBuilder,

    /*
    sim: FluidSim,
    c: DensitySim,
    m: DensitySim,
    y: DensitySim,
    k: DensitySim,
    */
    sim: CrystalSim,

    frame_count: usize,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let mut line_gb = GraphicsBuilder::new();
        let mut tri_gb = GraphicsBuilder::new();

        let mut rng = rand::thread_rng();

        let w = 200;
        let mut state = Array2D::new(w, w);
        state.data_mut().iter_mut().for_each(|b| *b = rng.gen_bool(0.5));

        let sim = CrystalSim::new(state);


        //draw_velocity_lines(&mut line_gb, sim.uv(), VELOCITY_Z);
        //draw_density(&mut tri_gb, c.density(), m.density(), y.density(), k.density(), DENSITY_Z);

        draw_crystal(&mut tri_gb, &sim.state, DENSITY_Z);

        //let line_verts = ctx.vertices(&line_gb.vertices, true)?;
        //let line_indices = ctx.indices(&line_gb.indices, true)?;

        let tri_verts = ctx.vertices(&tri_gb.vertices, true)?;
        let tri_indices = ctx.indices(&tri_gb.indices, true)?;

        let line_shader = ctx.shader(
            DEFAULT_VERTEX_SHADER,
            DEFAULT_FRAGMENT_SHADER,
            Primitive::Lines,
        )?;

        Ok(Self {
            /*
            line_verts,
            line_indices,
            line_gb,
            line_shader,
            */

            //c, m, y, k,

            tri_verts,
            tri_indices,
            tri_gb,

            sim,

            frame_count: 0,
        })
    }

    fn frame(&mut self, ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        // Modify
        self.frame_count += 1;

        self.sim.step(1000);

        // Draw
        //self.line_gb.clear();
        self.tri_gb.clear();

        draw_crystal(&mut self.tri_gb, &self.sim.state, DENSITY_Z);
        //draw_density(&mut self.tri_gb, self.c.density(), self.m.density(), self.y.density(), self.k.density(), DENSITY_Z);
        //draw_velocity_lines(&mut self.line_gb, self.sim.uv(), VELOCITY_Z);

        ctx.update_vertices(self.tri_verts, &self.tri_gb.vertices)?;
        //ctx.update_vertices(self.line_verts, &self.line_gb.vertices)?;

        // Render
        Ok(vec![
            DrawCmd::new(self.tri_verts).indices(self.tri_indices),
            /*DrawCmd::new(self.line_verts)
                .indices(self.line_indices)
                .shader(self.line_shader),*/
        ])
    }

    fn event(&mut self, ctx: &mut Context, platform: &mut Platform, event: Event) -> Result<()> {
        idek::simple_ortho_cam_ctx(ctx, platform);
        idek::close_when_asked(platform, &event);
        Ok(())
    }
}

struct CrystalSim {
    state: Array2D<bool>,
}

impl CrystalSim {
    pub fn new(state: Array2D<bool>) -> Self {
        Self { state }
    }

    pub fn step(&mut self, iters: usize) {
        let mut rng = rand::thread_rng();

        let neighbors = [
            (-1, 0),
            (1, 0),
            (0, -1),
            (0, 1),
        ];

        for _ in 0..iters {
            let x = rng.gen_range(1..self.state.width()-1);
            let y = rng.gen_range(1..self.state.height()-1);
            let (nx, ny) = neighbors.choose(&mut rng).unwrap();

            let a = (x, y);
            let b = (
                (x as i32 + nx) as usize,
                (y as i32 + ny) as usize,
            );

            let sa = self.state[a];
            let sb = self.state[b];
            if sa != sb {
                //self.state[a] = !sa;
                self.state[b] = !sb;
            }
        }

    }
}

fn draw_crystal(builder: &mut GraphicsBuilder, state: &Array2D<bool>, z: f32) {
    let cell_width = 2. / state.width() as f32;
    let cell_height = 2. / state.height() as f32;

    for i in 0..state.width() {
        let i_frac = (i as f32 / state.width() as f32) * 2. - 1.;
        for j in 0..state.height() {
            let j_frac = (j as f32 / state.height() as f32) * 2. - 1.;

            let c = state[(i, j)];
            let color = if c { [1.; 3] } else { [0.; 3] };

            let mut push = |dx: f32, dy: f32| {
                let pos = [i_frac + dx, j_frac + dy, z];
                builder.push_vertex(Vertex::new(pos, color))
            };

            let tl = push(0., 0.);
            let tr = push(cell_width, 0.);

            let bl = push(0., cell_height);
            let br = push(cell_width, cell_height);

            builder.push_indices(&[bl, tr, tl, bl, br, tr]);
        }
    }
}
