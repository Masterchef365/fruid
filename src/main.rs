use fruid::{Array2D, DensitySim, FluidSim};
use idek::{prelude::*, IndexBuffer};
mod graphics_builder;
use graphics_builder::GraphicsBuilder;
use rayon::prelude::*;

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

fn react(r: &mut f32, g: &mut f32, b: &mut f32, k: &mut f32, u: &mut f32, v: &mut f32) {
    //*u = (*k * 18.).sin() * *k * 8.;
    //*v = (*k * 19.).sin() * *k * 8.;

    //let d = *u * *u + *v * *v;

    let m = 0.2;
    let i = 2.0;
    let combust_rate = r.min(*g) / (1. + i*(*b - m).powf(2.)); // / (1. + *b);

    *r = -combust_rate;
    *g = -combust_rate;
    *b = combust_rate * 2.;
    *k = 0.;

    let kaboom = 10_000.;
    *u *= 1. + combust_rate * kaboom;
    *v *= 1. + combust_rate * kaboom;
}

fn reaction(
    c: &mut Array2D,
    m: &mut Array2D,
    y: &mut Array2D,
    k: &mut Array2D,
    u: &mut Array2D,
    v: &mut Array2D,
) {
    c.data_mut()
        .iter_mut()
        .zip(m.data_mut())
        .zip(y.data_mut())
        .zip(k.data_mut())
        .zip(u.data_mut())
        .zip(v.data_mut())
        .for_each(|(((((c, m), y), k), u), v)| react(c, m, y, k, u, v));
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let mut line_gb = GraphicsBuilder::new();
        let mut tri_gb = GraphicsBuilder::new();

        let mut sim = FluidSim::new(200, 200);

        let [mut c, mut m, mut y, mut k] =
            [(); 4].map(|_| DensitySim::new(sim.width(), sim.height()));

        let height = sim.height();
        let width = sim.width();
        let intensity = 15_000.;
        c.density_mut()[(width / 3, height / 2)] = intensity;
        m.density_mut()[(2 * width / 3, height / 2)] = intensity;
        //y.density_mut()[(4 * width / 5, height / 2)] = intensity;
        //k.density_mut()[(3 * width / 5, height / 2)] = intensity / 10.;

        sim.step(0.1, 0.0);
        c.step(sim.uv(), 0.1, 0.);
        m.step(sim.uv(), 0.1, 0.);
        y.step(sim.uv(), 0.1, 0.);
        k.step(sim.uv(), 0.1, 0.);

        dbg!(y.density()[(2 * width / 3, 2 * height / 3)]);

        draw_velocity_lines(&mut line_gb, sim.uv(), VELOCITY_Z);
        draw_density(
            &mut tri_gb,
            c.density(),
            m.density(),
            y.density(),
            k.density(),
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

            c,
            m,
            y,
            k,

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
        let time = self.frame_count as f32 / 12.; //ctx.start_time().elapsed().as_secs_f32();

        let (u, v) = self.sim.uv_mut();

        let center = (u.width() / 2, u.height() / 2);
        let x = center.0 as f32 * (-(time / 5.).cos() * 0.8 + 1.);
        u[(x as usize, center.1)] += time.cos() * 100.;

        // Step
        //self.c.density_mut().data_mut().fill(0.0);
        //self.m.density_mut().data_mut().fill(0.0);
        //self.y.density_mut().data_mut().fill(0.0);
        //self.k.density_mut().data_mut().fill(0.0);

        reaction(
            self.c.density_mut(),
            self.m.density_mut(),
            self.y.density_mut(),
            self.k.density_mut(),
            u,
            v,
        );

        let dt = 1e-2;
        let visc = 0.0;
        let diff = 0.0;

        self.sim.step(dt, visc);
        self.c.step(self.sim.uv(), dt, diff);
        self.m.step(self.sim.uv(), dt, diff);
        self.y.step(self.sim.uv(), dt, diff);
        self.k.step(self.sim.uv(), dt, diff);

        // Draw
        self.line_gb.clear();
        self.tri_gb.clear();

        draw_density(
            &mut self.tri_gb,
            self.c.density(),
            self.m.density(),
            self.y.density(),
            self.k.density(),
            DENSITY_Z,
        );
        draw_velocity_lines(&mut self.line_gb, self.sim.uv(), VELOCITY_Z);

        ctx.update_vertices(self.tri_verts, &self.tri_gb.vertices)?;
        ctx.update_vertices(self.line_verts, &self.line_gb.vertices)?;

        /*
        draw_density_image(
            &format!("{:05}.ppm", self.frame_count),
            self.c.density(),
            self.m.density(),
            self.y.density(),
            self.k.density(),
        );
        */

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

fn draw_density(
    builder: &mut GraphicsBuilder,
    c: &Array2D,
    m: &Array2D,
    y: &Array2D,
    k: &Array2D,
    z: f32,
) {
    let cell_width = 2. / c.width() as f32;
    let cell_height = 2. / c.height() as f32;

    for i in 0..c.width() {
        let i_frac = (i as f32 / c.width() as f32) * 2. - 1.;
        for j in 0..c.height() {
            let j_frac = (j as f32 / c.height() as f32) * 2. - 1.;

            // CMY dye
            let color = dye_colors(c[(i, j)], m[(i, j)], y[(i, j)], k[(i, j)]);

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

fn draw_density_image(path: &str, c: &Array2D, m: &Array2D, y: &Array2D, k: &Array2D) {
    let image: Vec<u8> = c
        .data()
        .iter()
        .zip(m.data())
        .zip(y.data())
        .zip(k.data())
        .map(|(((&c, &m), &y), &k)| dye_colors(c, m, y, k).map(|f| (f.clamp(0., 1.) * 256.) as u8))
        .flatten()
        .collect();
    write_ppm(path, &image, c.width()).expect("Failed to write image");
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

fn dye_colors(c: f32, m: f32, y: f32, k: f32) -> [f32; 3] {
    [c, m, y].map(|f| f + k)
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
