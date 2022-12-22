use std::{path::Path, io::{BufReader, BufWriter}, fs::File};

use anyhow::Result;
use fruid::FluidState;


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

    dbg!(input_image.width());
    dbg!(input_image.height());

    write_png("ok.png", &input_image)?;

    Ok(())
}

fn read_file(base: impl AsRef<Path>, number: usize) -> Result<FluidState> {
    let path = base.as_ref().join(format!("{}.bsl", number));
    let mut f = BufReader::new(File::open(path)?);
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

