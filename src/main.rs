use std::f32::consts::PI;

use fruid::FluidSim;
use idek::{prelude::*, IndexBuffer};
use idek_basics::Array2D;
use idek_basics::{idek, GraphicsBuilder};
use num_complex::Complex32;

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

    frame_count: usize,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let mut line_gb = GraphicsBuilder::new();
        let mut tri_gb = GraphicsBuilder::new();

        let mut sim = FluidSim::new(250, 250);
        //sim.step(0.1, 0.0, 10);

        draw_velocity_lines(&mut line_gb, sim.uv(), VELOCITY_Z);
        draw_density(&mut tri_gb, sim.smoke_mut(), DENSITY_Z);

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

            frame_count: 0,
        })
    }

    fn frame(&mut self, ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        // Modify
        self.frame_count += 1;
        let time = self.frame_count as f32 / 32.; //ctx.start_time().elapsed().as_secs_f32();

        let d = self.sim.smoke_mut();
        let center = (d.width() / 2, d.height() / 2);
        let x = center.0 as f32 * ((time / 5.).cos() + 1.) * 0.8;

        let (u, v) = self.sim.uv_mut();

        let pos = (x as usize, center.1);
        //let time = 3. * PI / 2.;
        u[pos] = -2500. * (time).cos();
        v[pos] = -2500. * (time).sin();

        let dt = 1e-2;
        let overstep = 1.9;

        let height = self.sim.height();
        let width = self.sim.width();
        self.sim.smoke_mut()[(2 * width / 5, height / 2)] += Complex32::from_polar(50., 0.);
        self.sim.smoke_mut()[(3 * width / 5, height / 2)] += Complex32::from_polar(50., PI / 2.);



        self.sim.step(dt, overstep, 15);

        // Draw
        self.line_gb.clear();
        self.tri_gb.clear();

        draw_density(&mut self.tri_gb, self.sim.smoke_mut(), DENSITY_Z);
        draw_velocity_lines(&mut self.line_gb, self.sim.uv(), VELOCITY_Z);

        ctx.update_vertices(self.tri_verts, &self.tri_gb.vertices)?;
        ctx.update_vertices(self.line_verts, &self.line_gb.vertices)?;

        // Render
        Ok(vec![
            DrawCmd::new(self.tri_verts).indices(self.tri_indices),
            DrawCmd::new(self.line_verts)
                .indices(self.line_indices)
                .shader(self.line_shader),
        ])
    }

    fn event(&mut self, ctx: &mut Context, platform: &mut Platform, event: Event) -> Result<()> {
        idek::simple_ortho_cam_ctx(ctx, platform);
        idek::close_when_asked(platform, &event);
        Ok(())
    }
}

fn draw_density(builder: &mut GraphicsBuilder, smoke: &Array2D<Complex32>, z: f32) {
    let cell_width = 2. / smoke.width() as f32;
    let cell_height = 2. / smoke.height() as f32;

    for i in 0..smoke.width() {
        let i_frac = (i as f32 / smoke.width() as f32) * 2. - 1.;
        for j in 0..smoke.height() {
            let j_frac = (j as f32 / smoke.height() as f32) * 2. - 1.;

            let k = smoke[(i, j)];

            let (r, theta) = k.to_polar();
            let theta = theta * 4.;
            let color = hsv_to_rgb(theta.to_degrees() + 180., 1., r);

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

fn draw_velocity_lines(b: &mut GraphicsBuilder, (u, v): (&Array2D<f32>, &Array2D<f32>), z: f32) {
    let cell_width = 2. / u.width() as f32;
    let cell_height = 2. / u.height() as f32;

    let step = 5;

    for i in (0..u.width()).step_by(step) {
        let i_frac = (i as f32 / u.width() as f32) * 2. - 1.;
        for j in (0..u.height()).step_by(step) {
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

/// https://gist.github.com/fairlight1337/4935ae72bcbcc1ba5c72
fn hsv_to_rgb(h: f32, s: f32, v: f32) -> [f32; 3] {
    let c = v * s; // Chroma
    let h_prime = (h / 60.0) % 6.0;
    let x = c * (1.0 - ((h_prime % 2.0) - 1.0).abs());
    let m = v - c;

    let (mut r, mut g, mut b);

    if 0. <= h_prime && h_prime < 1. {
        r = c;
        g = x;
        b = 0.0;
    } else if 1.0 <= h_prime && h_prime < 2.0 {
        r = x;
        g = c;
        b = 0.0;
    } else if 2.0 <= h_prime && h_prime < 3.0 {
        r = 0.0;
        g = c;
        b = x;
    } else if 3.0 <= h_prime && h_prime < 4.0 {
        r = 0.0;
        g = x;
        b = c;
    } else if 4.0 <= h_prime && h_prime < 5.0 {
        r = x;
        g = 0.0;
        b = c;
    } else if 5.0 <= h_prime && h_prime < 6.0 {
        r = c;
        g = 0.0;
        b = x;
    } else {
        r = 0.0;
        g = 0.0;
        b = 0.0;
    }

    r += m;
    g += m;
    b += m;

    [r, g, b]
}
