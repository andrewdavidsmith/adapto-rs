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
```
cargo build --release
```
in the root of the source directory, and then run
```
./target/release/adapto-rs
```
to see the command line arguments.
