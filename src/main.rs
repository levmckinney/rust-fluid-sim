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
    let mut fluid_sim = FluidSim::<ConstantForce>::new(SimConfig::default(), ConstantForce{ force: na::Vector3::new(0.0, -0.00098, 0.0)});
    fluid_sim.step();
}