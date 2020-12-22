# CSC 417 Final Project Lev Mckinney
## Setup instruction
This project is written in rust so to run it you will need the rust tool chain installed on your system.
Instruction for installing Rust can be found here https://www.rust-lang.org/tools/install.
If you are running ubuntu 20.04 you will also need to install some dependencies for bevy which 
you can is the game engine I use for rendering they can be installed with
```bash
$ sudo apt-get install pkg-config libx11-dev libasound2-dev libudev-dev
```
If you are using any other linux distribution see https://github.com/bevyengine/bevy/blob/master/docs/linux_dependencies.md.
Otherwise if you are running on OSX or Windows your good to go.
```bash
# make sure you have the rust tool chain on your path and run.
$ cargo run --release
```
This should build the code and start the simulation :D
## IDE for reading through the code
If you want to read the code with the ability to jump to definitions and get documentation when you hover
over things I recommend VSC with the `rust-lang.rust` and `matklad.rust-analyzer` extensions installed.
## Project structure
Physics simulation code in the project is held within `src/fluid_sim`
To start reading through start with `src/fluid_sim/mod.rs`