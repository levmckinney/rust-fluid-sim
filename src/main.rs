#[macro_use] extern crate itertools;

mod viewer;
mod fluid_sim;

use fluid_sim::{SimConfig};
use viewer::FluidSimViewer;

use bevy::{prelude::*};

fn main() {
    
    println!("To look around the simulation use the arrow keys.");
    println!("Press P view the pressures within the simulation grid.");
    println!("Press V to view the velocity components on the simulation grid.");
    App::build()
        .add_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(FluidSimViewer{ 
            config: SimConfig::default(),
        })
        .run();
}