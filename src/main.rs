use fruid::{Array2D, DensitySim, FluidSim};
use idek::{prelude::*, IndexBuffer};
mod graphics_builder;
use graphics_builder::GraphicsBuilder;

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
    r: DensitySim,
    g: DensitySim,
    b: DensitySim,

    frame_count: usize,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let mut line_gb = GraphicsBuilder::new();
        let mut tri_gb = GraphicsBuilder::new();

        let mut sim = FluidSim::new(250, 250);

        let [mut r, mut g, mut b] = [(); 3].map(|_| DensitySim::new(sim.width(), sim.height()));

        let height = sim.height();
        let width = sim.width();
        let intensity = 80. * (width * height) as f32;
        r.density_mut()[(width / 4, height / 2)] = intensity;
        g.density_mut()[(2 * width / 4, height / 2)] = intensity;
        b.density_mut()[(3 * width / 4, height / 2)] = intensity;

        sim.step(0.1, 0.0);
        r.step(sim.uv(), 0.1, 0.);
        g.step(sim.uv(), 0.1, 0.);
        b.step(sim.uv(), 0.1, 0.);

        draw_velocity_lines(&mut line_gb, sim.uv(), VELOCITY_Z);
        draw_density(
            &mut tri_gb,
            r.density(),
            g.density(),
            b.density(),
            DENSITY_Z,
        );

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

            r,
            g,
            b,

            tri_verts,
            tri_indices,
            tri_gb,

            sim,

            frame_count: 0,
        })
    }

    fn frame(&mut self, ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        // Modify
        let time = self.frame_count as f32 / 12.; //ctx.start_time().elapsed().as_secs_f32();

        let d = self.r.density_mut();
        let center = (d.width() / 2, d.height() / 2);
        let x = center.0 as f32 * ((time / 5.).cos() + 1.);

        let (u, v) = self.sim.uv_mut();

        let pos = (x as usize, center.1);
        u[pos] = -3500. * (time * 3.).cos();
        v[pos] = -3500. * (time * 3.).sin();

        // Step
        self.r.density_mut().data_mut().fill(0.0);
        self.g.density_mut().data_mut().fill(0.0);
        self.b.density_mut().data_mut().fill(0.0);

        let dt = 1e-2;
        let visc = 0.;
        let diff = 0.;

        self.sim.step(dt, visc);
        self.r.step(self.sim.uv(), dt, diff);
        self.g.step(self.sim.uv(), dt, diff);
        self.b.step(self.sim.uv(), dt, diff);

        // Draw
        self.line_gb.clear();
        self.tri_gb.clear();

        draw_density(
            &mut self.tri_gb,
            self.r.density(),
            self.g.density(),
            self.b.density(),
            DENSITY_Z,
        );

        draw_density_image(
            &format!("{}.ppm", self.frame_count),
            self.r.density(),
            self.g.density(),
            self.b.density(),
        );
        //draw_velocity_lines(&mut self.line_gb, self.sim.uv(), VELOCITY_Z);

        ctx.update_vertices(self.tri_verts, &self.tri_gb.vertices)?;
        //ctx.update_vertices(self.line_verts, &self.line_gb.vertices)?;

        self.frame_count += 1;

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

fn draw_density(builder: &mut GraphicsBuilder, r: &Array2D, g: &Array2D, b: &Array2D, z: f32) {
    let cell_width = 2. / r.width() as f32;
    let cell_height = 2. / r.height() as f32;

    for i in 0..r.width() {
        let i_frac = (i as f32 / r.width() as f32) * 2. - 1.;
        for j in 0..r.height() {
            let j_frac = (j as f32 / r.height() as f32) * 2. - 1.;

            // CMY dye
            let dye = [r[(i, j)], g[(i, j)], b[(i, j)]];

            let color = dye_colors(dye);

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

fn draw_density_image(path: &str, r: &Array2D, g: &Array2D, b: &Array2D) {
    let image: Vec<u8> = r
        .data()
        .iter()
        .zip(g.data())
        .zip(b.data())
        .map(|((&c, &m), &y)| dye_colors([c, m, y]).map(|f| (f.clamp(0., 1.) * 256.) as u8))
        .flatten()
        .collect();
    write_ppm(path, &image, r.width()).expect("Failed to write image");
}

pub fn write_ppm(path: &str, image: &[u8], width: usize) -> std::io::Result<()> {
    use std::io::Write;
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    let height = image.len() / (width * 3);
    assert_eq!(image.len() % 3, 0);
    assert_eq!(image.len() % (width * 3), 0);
    writer.write(format!("P6\n{} {}\n255\n", width, height).as_bytes())?;
    writer.write(image)?;
    Ok(())
}

fn dye_colors(dye: [f32; 3]) -> [f32; 3] {
    let total_dye = dye.iter().map(|d| d * d).sum::<f32>().sqrt();
    dye.map(|f| (1. - f) / (total_dye * 0.8).max(1.))
        .map(|d| (d + 2.).log2())
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
