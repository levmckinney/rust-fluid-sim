mod viewer;
use bevy::{prelude::*};
use nalgebra as na;
use viewer::FluidSimViewer;

fn main() {
    App::build()
        .add_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(FluidSimViewer{ 
            sim_volume: viewer::SimVolume {
                num_vis_particles: 32,
                center:na::Vector3::new(0.0, 1.0, 0.0), 
                size: na::Vector3::new(1.0, 1.0, 1.0)
            }})
        .run();
}