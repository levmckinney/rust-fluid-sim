#[macro_use] extern crate itertools;

mod viewer;
mod fluid_sim;

use fluid_sim::{FluidSim, SimConfig};
use fluid_sim::forces::ConstantForce;
use viewer::FluidSimViewer;

use bevy::{prelude::*};
use nalgebra as na;

fn main() {
    App::build()
        .add_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(FluidSimViewer{ 
            config: SimConfig::default(),
        })
        .run();
}