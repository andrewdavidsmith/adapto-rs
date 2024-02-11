# adapto-rs

Remove adaptors from short-read sequencing data. I needed a faster
tool to remove adaptor sequences from reads for my own in my own data
analysis. I wrote this in Rust because I wanted to experiment with the
libraries for parallelism. It turns out this application does not need
much in the way of parallelism. I think this code is reasonably fast.
Unlike the tools I was using previously, `adapto-rs` does not use more
cores concurrently than the number of threads specified by the user,
which is an absolute requirement of my computing environment.

If you have `cargo` installed, you can build this code by doing:
```console
cargo build --release
```
in the root of the source directory, and then run
```console
./target/release/adapto-rs
```
to see the command line arguments.

If you need a static linked binary, which I currently do for a cluster
(x86, linux) environment that I do not manage, you can build one on
your own Ubuntu machine by doing the following:
```console
sudo apt install musl-tools
rustup target add x86_64-unknown-linux-musl
RUSTFLAGS='-C target-feature=+crt-static' cargo build --release --target x86_64-unknown-linux-musl
```
I think I got a smaller and possibly faster binary with this instead:
```console
RUSTFLAGS='-C target-feature=+crt-static -C codegen-units=1 -C lto -C embed-bitcode=yes' cargo build --release --target x86_64-unknown-linux-musl
```
You will then find the binary in the
`./target/x86_64-unknown-linux-musl/release/` subdirectory.
