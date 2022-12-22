cargo run --release --bin fruid -- --output $1
cargo run --release --bin animate -- --record $1 --image $2
ffmpeg -i $1/%4d.png -framerate 60 $1/out.mp4
vlc $1/out.mp4
