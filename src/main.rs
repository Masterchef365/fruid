use std::{collections::HashSet, f32::consts::PI};

use fruid::{FluidSim, SmokeSim};
use idek::{prelude::*, IndexBuffer};
use idek_basics::{
    idek::{self, nalgebra::SimdBool},
    Array2D, GraphicsBuilder,
};

fn main() -> Result<()> {
    //plot_fill_circle((0,0), 0, |pt| {dbg!(pt);});
    //Ok(())

    launch::<_, TriangleApp>(Settings::default().vr_if_any_args())
}

type Color = [f32; 3];

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
    life: ParticleLife,

    frame_count: usize,
}

struct ParticleLife {
    smoke: Vec<SmokeSim>,
    behaviours: Array2D<Behaviour>,
    colors: Vec<Color>,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        // Set up fluid sim
        let w = 150;
        let sim = FluidSim::new(w, w);

        let height = sim.height();
        let width = sim.width();

        // Decide behaviours and colors
        let n = 2;
        let colors: Vec<Color> = (0..n)
            .map(|_| hsv_to_rgb(rand::random::<f32>() * 360., 1., 1.))
            .collect();
        let behaviours: Vec<Behaviour> = (0..n * n)
            .map(|_| {
                Behaviour::default().with_inter_strength((rand::random::<f32>() * 2. - 1.) * 5.)
            })
            .collect();
        let behaviours = Array2D::from_array(n, behaviours);

        // Create life!
        let mut life = ParticleLife::new(w, w, behaviours, colors);

        // Place a dot of smoke
        for smoke in life.smoke_mut() {
            for _ in 0..100 {
                let intensity = 10.;

                let area_width = width - 4;
                let area_height = height - 4;

                let x = area_width * (rand::random::<usize>() % area_width) / area_width + 2;
                let y = area_height * (rand::random::<usize>() % area_height) / area_height + 2;
                smoke.smoke_mut()[(x, y)] = intensity;
            }

            /*
                    let radius = 50;
                    let mut ctr = 0;
                    plot_fill_circle(, radius, |pt| {
                        if let Some(coord) = box_coord((w, w), (x, y)) {
            let coord =                (w as i32 / 2, w as i32 / 2);
                            smoke.smoke_mut()[coord] = 1.;
                        }
                    });
                    */
        }

        // Set up line buffer
        let mut line_gb = GraphicsBuilder::new();
        draw_velocity_lines(&mut line_gb, sim.uv(), VELOCITY_Z);

        let line_verts = ctx.vertices(&line_gb.vertices, true)?;
        let line_indices = ctx.indices(&line_gb.indices, true)?;

        // Set up triangle buffer
        let mut tri_gb = GraphicsBuilder::new();
        draw_density(&mut tri_gb, &life.smoke, &life.colors, DENSITY_Z);

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
            life,

            frame_count: 0,
        })
    }

    fn frame(&mut self, ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        /*
        // Modify
        self.frame_count += 1;
        let time = self.frame_count as f32 / 120.; //ctx.start_time().elapsed().as_secs_f32();

        let d = self.life.smoke_mut()[0].smoke_mut();
        let center = (d.width() / 2, d.height() / 2);
        let x = center.0 as f32 * ((time / 5.).cos() + 1.) * 0.8;

        let (u, v) = self.sim.uv_mut();

        let pos = (x as usize, center.1);
        let pos = center;
        //let time = 3. * PI / 2.;
        u[pos] = -450. * (time).cos();
        v[pos] = -450. * (time).sin();
        */

        let dt = 1e-2;
        let overstep = 1.9;

        particle_life(
            &self.life.behaviours,
            &self.life.smoke,
            self.sim.uv_mut(),
            1.,
        );

        self.sim.step(dt, overstep, 15);
        for smoke in self.life.smoke_mut() {
            smoke.advect(self.sim.uv(), dt);
        }

        // Draw
        self.line_gb.clear();
        self.tri_gb.clear();

        draw_density(
            &mut self.tri_gb,
            &mut self.life.smoke,
            &mut self.life.colors,
            DENSITY_Z,
        );
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

impl ParticleLife {
    pub fn new(
        width: usize,
        height: usize,
        behaviours: Array2D<Behaviour>,
        colors: Vec<Color>,
    ) -> Self {
        let smoke = (0..colors.len())
            .map(|_| SmokeSim::new(width, height))
            .collect();

        Self {
            smoke,
            behaviours,
            colors,
        }
    }

    pub fn smoke_mut(&mut self) -> &mut [SmokeSim] {
        &mut self.smoke
    }
}

fn draw_density(builder: &mut GraphicsBuilder, smoke: &[SmokeSim], colors: &[Color], z: f32) {
    let width = smoke[0].smoke().width();
    let height = smoke[0].smoke().height();

    let cell_width = 2. / width as f32;
    let cell_height = 2. / height as f32;

    for i in 0..width {
        let i_frac = (i as f32 / width as f32) * 2. - 1.;
        for j in 0..height {
            let j_frac = (j as f32 / height as f32) * 2. - 1.;

            // Sum colors
            let color = smoke
                .iter()
                .zip(colors)
                .map(|(smoke, color)| color.into_iter().map(|&c| c * smoke.smoke()[(i, j)]))
                .fold([0.; 3], |mut acc, x| {
                    acc.iter_mut().zip(x).for_each(|(acc, x)| *acc += x);
                    acc
                });

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

#[derive(Clone, Copy, Debug)]
pub struct Behaviour {
    /// Magnitude of the default repulsion force
    pub default_repulse: f32,
    /// Zero point between default repulsion and particle interaction (0 to 1)
    pub inter_threshold: f32,
    /// Interaction peak strength
    pub inter_strength: f32,
    /// Maximum distance of particle interaction (0 to 1)
    pub inter_max_dist: f32,
}

impl Behaviour {
    /// Returns the force on this particle
    fn interact(&self, dist: f32) -> f32 {
        if dist < self.inter_threshold {
            let f = dist / self.inter_threshold;
            (1. - f) * -self.default_repulse
        } else if dist > self.inter_max_dist {
            0.0
        } else {
            let x = dist - self.inter_threshold;
            let x = x / (self.inter_max_dist - self.inter_threshold);
            let x = x * 2. - 1.;
            let x = 1. - x.abs();
            x * self.inter_strength
        }
    }
}

/// https://gist.github.com/fairlight1337/4935ae72bcbcc1ba5c72
fn hsv_to_rgb(h: f32, s: f32, v: f32) -> Color {
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

impl Default for Behaviour {
    fn default() -> Self {
        Self {
            default_repulse: 1.,
            inter_threshold: 3.,
            inter_strength: 1.,
            inter_max_dist: 5.0,
        }
    }
}

impl Behaviour {
    pub fn with_inter_strength(mut self, inter_strength: f32) -> Self {
        self.inter_strength = inter_strength;
        self
    }
}

fn particle_life(
    behaviours: &Array2D<Behaviour>,
    smokes: &[SmokeSim],
    (u, v): (&mut Array2D<f32>, &mut Array2D<f32>),
    dt: f32,
) {
    let w = u.width();
    let h = u.height();

    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let (dx, dy) = calc_force((x, y), behaviours, smokes);
            u[(x, y)] += dx * dt;
            v[(x, y)] += dy * dt;
        }
    }
}

fn calc_force(
    center: Coord<usize>,
    behaviours: &Array2D<Behaviour>,
    smokes: &[SmokeSim],
) -> (f32, f32) {
    let width = smokes[0].smoke().width();
    let height = smokes[0].smoke().height();

    let mut force_x = 0.;
    let mut force_y = 0.;

    /*
    let circles: Vec<Vec<Coord<i32>>> = behaviours
        .data()
        .iter()
        .map(|b| plot_fill_circle(b.inter_max_dist as i32).collect())
        .collect();
    let circles = Array2D::from_array(width, circles);
    */

    // For each possible behaviour (interaction)...
    for i in 0..smokes.len() {
        let i_smoke = &smokes[i];
        for j in 0..smokes.len() {
            let j_smoke = &smokes[j];

            let behav = behaviours[(i, j)];
            let radius = behav.inter_max_dist as i32;

            // For each point in the circle around this point
            for (dx, dy) in plot_fill_circle(radius) {
                // Skip thyself
                if (dx, dy) == (0, 0) {
                    continue;
                }

                let (cx, cy) = center;
                let pt = (cx as i32 + dx, cy as i32 + dy);

                if let Some(coord) = box_coord((width, height), pt) {
                    // Calculate relative distance
                    let dist_sq = dx * dx + dy * dy;
                    let dist = (dist_sq as f32).sqrt();

                    // Calculate force magnitude
                    let total_smoke = i_smoke.smoke()[center] * j_smoke.smoke()[coord];
                    let inter = behav.interact(dist);
                    let mag = inter * total_smoke;

                    // Calculate forces
                    force_x += mag * dx as f32 / dist.max(1.);
                    force_y += mag * dy as f32 / dist.max(1.);
                }
            }
        }
    }

    (force_x, force_y)
}

type Coord<T> = (T, T);

fn box_coord((width, height): Coord<usize>, (x, y): Coord<i32>) -> Option<Coord<usize>> {
    let f = |v: i32, len: usize| (v >= 0).then(|| v as usize).filter(|&v| v < len);
    f(x, width).zip(f(y, height))
}

fn plot_circle((x0, y0): Coord<i32>, radius: i32, mut plot: impl FnMut(Coord<i32>)) {
    let mut x = radius;
    let mut y = 0;
    let mut err = 0;

    while x >= y {
        // Draw the top and bottom points of the circle.
        plot((x0 + x, y0 + y));
        plot((x0 - x, y0 + y));
        plot((x0 + x, y0 - y));
        plot((x0 - x, y0 - y));

        // Draw a horizontal line between the top and bottom points.
        if x != y {
            plot((x0 + y, y0 + x));
            plot((x0 - y, y0 + x));
            plot((x0 + y, y0 - x));
            plot((x0 - y, y0 - x));
        }

        // Calculate the next point on the circumference.
        y += 1;
        err += 1 + 2 * y;
        if 2 * (err - x) + 1 > 0 {
            x -= 1;
            err += 1 - 2 * x;
        }
    }
}

fn plot_fill_circle(radius: i32) -> impl Iterator<Item = Coord<i32>> {
    // TODO: Use midpoint circle algo smh
    (-radius..=radius)
        .map(move |y| {
            (-radius..=radius)
                .filter_map(move |x| (x * x + y * y < radius * radius).then(|| (x, y)))
        })
        .flatten()
}
