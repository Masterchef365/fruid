use std::{path::Path, io::{BufReader, BufWriter}, fs::File};

use anyhow::Result;
use fruid::{FluidState, interp, advect};


use std::path::PathBuf;
use structopt::StructOpt;

type Image = fruid::array2d::Array2D<[u8; 3]>;

#[derive(Debug, StructOpt, Default)]
#[structopt(name = "example", about = "An example of StructOpt usage.")]
struct Opt {
    #[structopt(short, long)]
    record: PathBuf,

    #[structopt(short, long)]
    image: PathBuf,
}

fn main() -> Result<()> {
    let args = Opt::from_args();
    let n = last_number(&args.record);

    let input_image = load_png(&args.image)?;
    let example = read_file(&args.record, 0)?;

    let mut parts = init_particles(&example, &input_image);

    let dt = 0.1;

    for i in 0..n {
        println!("{}/{}", i+1, n);
        let path = args.record.join(format!("{}.png", n - i));
        let state = read_file(&args.record, n - i)?;
        step_particles(&state, &mut parts, dt);
        let out_img = particle_image(&example, &parts, &input_image);
        write_png(&path, &out_img)?;
    }

    Ok(())
}

fn particle_image(example: &FluidState, parts: &[[f32; 2]], input_img: &Image) -> Image {
    let mut output_img = Image::new(input_img.width(), input_img.height());

    let w = example.u.width() as f32 - 1.;
    let h = example.v.height() as f32 - 1.;

    let mut n = 0;
    for y in 0..input_img.height() {
        for x in 0..input_img.width() {
            let pixel = input_img[(x, y)];

            let [px, py] = parts[n];

            let ux = ((px - 0.5) / w) * input_img.width() as f32;
            let uy = ((py - 0.5) / h) * input_img.height() as f32;

            let x_bound = ux.floor() > 0. && ux.floor() < input_img.width() as f32;
            let y_bound = uy.floor() > 0. && uy.floor() < input_img.height() as f32;

            if x_bound && y_bound {
                let pos = (ux as usize, uy as usize);
                output_img[pos] = pixel;
            }

            n += 1;
        }
    }

    output_img
}

fn step_particles(state: &FluidState, parts: &mut [[f32; 2]], dt: f32) {
    let w = state.u.width() as f32 - 1.;
    let h = state.v.height() as f32 - 1.;

    for [x, y] in parts {
        if x.floor() > 0.5 && x.floor() < w - 1. && y.floor() > 0.5 && y.floor() < h - 1. {
            let (ax, ay) = advect(&state.u, &state.v, *x, *y, dt);
            *x = ax;
            *y = ay;
        }
    }
}

fn init_particles(example: &FluidState, image: &Image) -> Vec<[f32; 2]> {
    let w = example.u.width() as f32 - 1.;
    let h = example.v.height() as f32 - 1.;

    let mut parts = vec![];

    for y in 0..image.height() {
        for x in 0..image.width() {
            let px = (x as f32 / image.width() as f32) * w + 0.5;
            let py = (y as f32 / image.height() as f32) * h + 0.5;
            parts.push([px, py]);
        }
    }
    
    parts
}

fn read_file(base: impl AsRef<Path>, number: usize) -> Result<FluidState> {
    let path = base.as_ref().join(format!("{}.bsl", number));
    let f = BufReader::new(File::open(path)?);
    Ok(bincode::deserialize_from(f)?)
}

fn last_number(base: impl AsRef<Path>) -> usize {
    let mut n = 0;
    loop {
        let p = base.as_ref().join(format!("{}.bsl", n));
        if !p.is_file() {
            return n - 1;
        }
        n += 1;
    }
}

fn load_png(path: impl AsRef<Path>) -> Result<Image> {
    let decoder = png::Decoder::new(File::open(path)?);

    let mut reader = decoder.read_info()?;
    let mut buf = vec![0; reader.output_buffer_size()];

    let info = reader.next_frame(&mut buf)?;

    assert!(info.color_type == png::ColorType::Rgb);
    assert!(info.bit_depth == png::BitDepth::Eight);

    buf.truncate(info.buffer_size());

    let mut data = vec![];
    for pixel in buf.chunks_exact(3) {
        let mut k = [0; 3];
        k.copy_from_slice(pixel);
        data.push(k);
    }

    Ok(Image::from_array(info.width as _, data))
}

fn write_png(path: impl AsRef<Path>, image: &Image) -> Result<()> {
    let file = File::create(path)?;
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, image.width() as _, image.height() as _); // Width is 2 pixels and height is 1.
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);

    let mut writer = encoder.write_header()?;
    let data: Vec<u8> = image.data().iter().copied().flatten().collect();
    writer.write_image_data(&data)?; // Save

    Ok(())
}

