use idek_basics::{idek, GraphicsBuilder};
use fruid::{Array2D, FluidSim, DensitySim};
use idek::{prelude::*, IndexBuffer};

fn main() -> Result<()> {
    launch::<_, TriangleApp>(Settings::default().vr_if_any_args())
}

const DENSITY_Z: f32 = 0.5;
const VELOCITY_Z: f32 = 0.0;

struct TriangleApp {
    line_verts: VertexBuffer,
    line_indices: IndexBuffer,
    line_gb: GraphicsBuilder,
    line_shader: Shader,

    tri_verts: VertexBuffer,
    tri_indices: IndexBuffer,
    tri_gb: GraphicsBuilder,

    sim: FluidSim,
    c: DensitySim,
    m: DensitySim,
    y: DensitySim,
    k: DensitySim,

    frame_count: usize,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let mut line_gb = GraphicsBuilder::new();
        let mut tri_gb = GraphicsBuilder::new();

        let mut sim = FluidSim::new(250, 250);

        let [mut c, mut m, mut y, mut k] = [(); 4].map(|_| DensitySim::new(sim.width(), sim.height()));

        let height = sim.height();
        let width = sim.width();
        let intensity = 80. * (width * height) as f32;
        c.density_mut()[(width / 5, height / 2)] = intensity;
        k.density_mut()[(2 * width / 5, height / 2)] = intensity / 100.;
        m.density_mut()[(3 * width / 5, height / 2)] = intensity;
        y.density_mut()[(4 * width / 5, height / 2)] = intensity;

        sim.step(0.1, 0.0);
        c.step(sim.uv(), 0.1, 0.);
        m.step(sim.uv(), 0.1, 0.);
        y.step(sim.uv(), 0.1, 0.);
        k.step(sim.uv(), 0.1, 0.);

        dbg!(y.density()[(2 * width / 3, 2 * height / 3)]);

        draw_velocity_lines(&mut line_gb, sim.uv(), VELOCITY_Z);
        draw_density(&mut tri_gb, c.density(), m.density(), y.density(), k.density(), DENSITY_Z);

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

            c, m, y, k,

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
        let time = self.frame_count as f32 / 12.;//ctx.start_time().elapsed().as_secs_f32();

        let d = self.c.density_mut();
        let center = (d.width() / 2, d.height() / 2);
        let x = center.0 as f32 * ((time / 5.).cos() + 1.);

        let (u, v) = self.sim.uv_mut();

        let pos = (x as usize, center.1);
        u[pos] = -4500. * (time * 3.).cos();
        v[pos] = -4500. * (time * 3.).sin();

        // Step
        self.c.density_mut().data_mut().fill(0.0);
        self.m.density_mut().data_mut().fill(0.0);
        self.y.density_mut().data_mut().fill(0.0);
        self.k.density_mut().data_mut().fill(0.0);

        let dt = 1e-2;
        let visc = 0.;
        let diff = 0.;

        self.sim.step(dt, visc);
        self.c.step(self.sim.uv(), dt, diff);
        self.m.step(self.sim.uv(), dt, diff);
        self.y.step(self.sim.uv(), dt, diff);
        self.k.step(self.sim.uv(), dt, diff);

        // Draw
        self.line_gb.clear();
        self.tri_gb.clear();

        draw_density(&mut self.tri_gb, self.c.density(), self.m.density(), self.y.density(), self.k.density(), DENSITY_Z);
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

fn draw_density(builder: &mut GraphicsBuilder, c: &Array2D, m: &Array2D, y: &Array2D, k: &Array2D, z: f32) {
    let cell_width = 2. / c.width() as f32;
    let cell_height = 2. / c.height() as f32;

    for i in 0..c.width() {
        let i_frac = (i as f32 / c.width() as f32) * 2. - 1.;
        for j in 0..c.height() {
            let j_frac = (j as f32 / c.height() as f32) * 2. - 1.;

            // CMY dye
            let k = k[(i, j)];
            let cmy = [
                c[(i, j)],
                m[(i, j)],
                y[(i, j)],
            ];

            let total_dye = cmy.iter().map(|d| d * d).sum::<f32>().sqrt();

            let color = cmy
                .map(|f| (1. - f - k) / total_dye.max(1.))
                .map(|d| (d + 2.).log2());

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

fn draw_velocity_lines(b: &mut GraphicsBuilder, (u, v): (&Array2D, &Array2D), z: f32) {
    let cell_width = 2. / u.width() as f32;
    let cell_height = 2. / u.height() as f32;

    for i in 0..u.width() {
        let i_frac = (i as f32 / u.width() as f32) * 2. - 1.;
        for j in 0..u.height() {
            let j_frac = (j as f32 / u.height() as f32) * 2. - 1.;

            let u = u[(i, j)];
            let v = v[(i, j)];

            let speed = (u.powf(2.) + v.powf(2.)).sqrt();

            let color = [speed; 3];

            let mut push = |x: f32, y: f32| {
                let pos = [x, y, z];
                b.push_vertex(Vertex::new(pos, color))
            };

            let tail_x = i_frac + cell_width / 2.;
            let tail_y = j_frac + cell_height / 2.;
            let tail = push(tail_x, tail_y);

            let len = cell_height * 2. / speed;
            let tip = push(tail_x + u * len, tail_y + v * len);

            b.push_indices(&[tip, tail]);
        }
    }
}
