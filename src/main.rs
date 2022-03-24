use fruid::{Array3D, DensitySim, FluidSim};
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

    point_verts: VertexBuffer,
    point_indices: IndexBuffer,
    point_gb: GraphicsBuilder,
    point_shader: Shader,

    sim: FluidSim,
    c: DensitySim,
    m: DensitySim,
    y: DensitySim,
    black: DensitySim,

    frame_count: usize,

    camera: MultiPlatformCamera,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, platform: &mut Platform, _: ()) -> Result<Self> {
        let mut line_gb = GraphicsBuilder::new();
        let mut point_gb = GraphicsBuilder::new();

        let dim = 65;
        let mut sim = FluidSim::new(dim, dim, dim);

        let [mut c, mut m, mut y, mut black] =
            [(); 4].map(|_| DensitySim::new(sim.width(), sim.height(), sim.length()));

        let height = sim.height();
        let width = sim.width();
        let length = sim.length();

        let intensity = 1. * (width * height * length) as f32;
        //c.density_mut()[(width / 5, height / 2, length / 2)] = intensity;
        black.density_mut()[(2 * width / 5, height / 2, length / 2)] = intensity;
        //m.density_mut()[(3 * width / 5, height / 2, length / 2)] = intensity;
        //y.density_mut()[(4 * width / 5, height / 2, length / 2)] = intensity;

        sim.step(0.1, 0.0);
        c.step(sim.uvw(), 0.1, 0.);
        m.step(sim.uvw(), 0.1, 0.);
        y.step(sim.uvw(), 0.1, 0.);
        black.step(sim.uvw(), 0.1, 0.);

        dbg!(y.density()[(2 * width / 3, 2 * width / 3, 2 * height / 3)]);

        draw_velocity_lines(&mut line_gb, sim.uvw());
        draw_density(
            &mut point_gb,
            c.density(),
            m.density(),
            y.density(),
            black.density(),
        );

        let line_verts = ctx.vertices(&line_gb.vertices, true)?;
        let line_indices = ctx.indices(&line_gb.indices, true)?;

        let point_verts = ctx.vertices(&point_gb.vertices, true)?;
        let point_indices = ctx.indices(&point_gb.indices, true)?;

        let line_shader = ctx.shader(
            DEFAULT_VERTEX_SHADER,
            DEFAULT_FRAGMENT_SHADER,
            Primitive::Lines,
            Blend::Opaque,
        )?;

        let point_shader = ctx.shader(
            include_bytes!("shaders/unlit.vert.spv"),
            DEFAULT_FRAGMENT_SHADER,
            Primitive::Points,
            Blend::Additive,
        )?;

        Ok(Self {
            camera: MultiPlatformCamera::new(platform),

            line_verts,
            line_indices,
            line_gb,
            line_shader,

            c,
            m,
            y,
            black,

            point_verts,
            point_indices,
            point_gb,
            point_shader,

            sim,

            frame_count: 0,
        })
    }

    fn frame(&mut self, ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        // Modify
        self.frame_count += 1;
        let time = self.frame_count as f32 / 12.; //ctx.start_time().elapsed().as_secs_f32();

        let d = self.c.density_mut();
        let center = (d.width() / 2, d.height() / 2, d.length() / 2);
        let x = center.0 as f32 * ((time / 5.).cos() + 1.);

        let (u, v, w) = self.sim.uvw_mut();

        let pos = (x as usize, center.1, center.2);
        u[pos] = -4500. * (time * 3.).cos();
        v[pos] = -4500. * (time * 3.).sin();
        w[pos] = -4500. * (time * 2.).sin();

        // Step
        self.c.density_mut().data_mut().fill(0.0);
        self.m.density_mut().data_mut().fill(0.0);
        self.y.density_mut().data_mut().fill(0.0);
        self.black.density_mut().data_mut().fill(0.0);

        let dt = 1e-2;
        let visc = 0.0001;
        let diff = 0.;

        self.sim.step(dt, visc);
        //self.c.step(self.sim.uvw(), dt, diff);
        //self.m.step(self.sim.uvw(), dt, diff);
        //self.y.step(self.sim.uvw(), dt, diff);
        self.black.step(self.sim.uvw(), dt, diff);

        // Draw
        self.line_gb.clear();
        self.point_gb.clear();

        draw_density(&mut self.point_gb, self.c.density(), self.m.density(), self.y.density(), self.black.density());
        draw_velocity_lines(&mut self.line_gb, self.sim.uvw());

        ctx.update_vertices(self.point_verts, &self.point_gb.vertices)?;
        ctx.update_vertices(self.line_verts, &self.line_gb.vertices)?;

        // Render
        Ok(vec![
            DrawCmd::new(self.point_verts)
                .indices(self.point_indices)
                .shader(self.point_shader),
            /*DrawCmd::new(self.line_verts)
                .indices(self.line_indices)
                .shader(self.line_shader)*/
        ])
    }

    fn event(
        &mut self,
        ctx: &mut Context,
        platform: &mut Platform,
        mut event: Event,
    ) -> Result<()> {
        if self.camera.handle_event(&mut event) {
            ctx.set_camera_prefix(self.camera.get_prefix())
        }
        idek::close_when_asked(platform, &event);
        Ok(())
    }
}

fn draw_density(
    builder: &mut GraphicsBuilder,
    c: &Array3D,
    m: &Array3D,
    y: &Array3D,
    black: &Array3D,
) {
    for i in 0..c.width() {
        let i_frac = (i as f32 / c.width() as f32) * 2. - 1.;
        for j in 0..c.height() {
            let j_frac = (j as f32 / c.height() as f32) * 2. - 1.;
            for k in 0..c.length() {
                let k_frac = (k as f32 / c.length() as f32) * 2. - 1.;

                // CMY dye
                let black = black[(i, j, k)];
                let cmy = [c[(i, j, k)], m[(i, j, k)], y[(i, j, k)]];

                //let total_dye = cmy.iter().map(|d| d * d).sum::<f32>().sqrt();

                let color = cmy.map(|f| f + black);
                //.map(|d| (d + 2.).log2());

                let idx = builder.push_vertex(Vertex::new([i_frac, j_frac, k_frac], color));
                builder.push_indices(&[idx]);
            }
        }
    }
}

fn draw_velocity_lines(b: &mut GraphicsBuilder, (u, v, w): (&Array3D, &Array3D, &Array3D)) {
    let cell_width = 2. / u.width() as f32;
    let cell_height = 2. / u.height() as f32;
    let cell_length = 2. / u.length() as f32;

    for i in 0..u.width() {
        let i_frac = (i as f32 / u.width() as f32) * 2. - 1.;
        for j in 0..u.height() {
            let j_frac = (j as f32 / u.height() as f32) * 2. - 1.;
            for k in 0..u.length() {
                let k_frac = (k as f32 / u.length() as f32) * 2. - 1.;

                let u = u[(i, j, k)];
                let v = v[(i, j, k)];
                let w = w[(i, j, k)];

                let speed = (u.powf(2.) + v.powf(2.)).sqrt();

                let color = [speed; 3];

                let mut push = |x: f32, y: f32, z: f32| {
                    let pos = [x, y, z];
                    b.push_vertex(Vertex::new(pos, color))
                };

                let tail_x = i_frac + cell_width / 2.;
                let tail_y = j_frac + cell_height / 2.;
                let tail_z = k_frac + cell_length / 2.;
                let tail = push(tail_x, tail_y, tail_z);

                let len = cell_height * 1. / speed;
                let tip = push(tail_x + u * len, tail_y + v * len, tail_z + w * len);

                b.push_indices(&[tip, tail]);
            }
        }
    }
}
