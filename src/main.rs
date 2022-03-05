use idek::{prelude::*, IndexBuffer};

fn main() -> Result<()> {
    launch::<_, TriangleApp>(Settings::default().vr_if_any_args())
}

const QUAD_VERTS: [Vertex; 4] = [
    Vertex {
        pos: [-1., -1., 0.],
        color: [0., 0., 0.],
    },
    Vertex {
        pos: [1., -1., 0.],
        color: [1., 0., 0.],
    },
    Vertex {
        pos: [1., 1., 0.],
        color: [1., 1., 0.],
    },
    Vertex {
        pos: [-1., 1., 0.],
        color: [0., 1., 0.],
    },
];

const QUAD_INDICES: [u32; 12] = [
    // Facing toward the camera
    3, 1, 0, 3, 2, 1, // Facing away
    0, 1, 3, 1, 2, 3,
];

struct TriangleApp {
    verts: VertexBuffer,
    indices: IndexBuffer,
}

impl App for TriangleApp {
    fn init(ctx: &mut Context, _: &mut Platform, _: ()) -> Result<Self> {
        let verts = ctx.vertices(&QUAD_VERTS, false)?;
        let indices = ctx.indices(&QUAD_INDICES, false)?;

        Ok(Self {
            verts,
            indices,
        })
    }

    fn frame(&mut self, _ctx: &mut Context, _: &mut Platform) -> Result<Vec<DrawCmd>> {
        Ok(vec![DrawCmd::new(self.verts)
            .indices(self.indices)])
    }

    fn event(
        &mut self,
        ctx: &mut Context,
        platform: &mut Platform,
        event: Event,
    ) -> Result<()> {
        idek::simple_ortho_cam_ctx(ctx, platform);
        idek::close_when_asked(platform, &event);
        Ok(())
    }
}

