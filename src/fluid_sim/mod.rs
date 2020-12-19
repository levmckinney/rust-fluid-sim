pub mod particle;
pub mod mac_grid;
pub mod forces;
//pub(crate) mod conjugate_gradient;
pub mod utils;

use particle::{Particle, Particles};
use forces::Force;
use mac_grid::MACGrid;
use utils::Grid;
use nalgebra as na;

/// density - This is the number of cells along each axis. 
/// dims - The dimension of the sim volume (depth, height, width).
/// Note that all sim volumes start from 0, 0, 0.
#[derive(Clone, Copy)]
pub struct SimConfig {
    delta_time: f32,
    density: f32,
    grid: Grid
}

impl Default for SimConfig {
    fn default() -> Self {
        SimConfig {
            delta_time: 1.0/60.0,
            density: 1.225, // density of air in Kg per meter cubed
            grid: Grid {
                dims: na::Vector3::new(1.0, 1.0, 1.0),
                grid_shape: (2, 2, 2)
            }
        }
    }
}

pub struct FluidSim<T> where T: Force {
    dt: f32,
    mac_grid: MACGrid,
    pub particles: Particles,
    force: T
}

/// Fluid simulator
impl<T> FluidSim <T> where T: Force{

    /// Arguments:
    /// * config - The configuration for the simulation
    /// Returns:
    /// * A newly constructed fluid sim
    pub fn new(config: SimConfig, force: T) -> Self {
        Self {
            dt: config.delta_time,
            mac_grid: MACGrid::new(&config.grid, config.density),
            particles: Particles::filled(&config.grid),
            force: force
        }
    }

    /// Arguments:
    /// * position - the position to sample at
    /// Returns:
    /// * the velocity at said position if it is within the simvolume otherwise
    /// * it returns nothing.
    pub fn get_velocity(pos: na::Vector3<f32>) -> Option<na::Vector3<f32>> {
        unimplemented!()
    }

    /// Arguments:
    /// * dt - amount of time to step the fluid simulation
    /// Results:
    /// * updates the state of the simulation by one time step of length dt.
    pub fn step(&mut self) {
        // TODO: Transfer particle velocities to grid
        self.particles.advect(self.dt);
        self.particles.apply_force(&self.force, self.dt)
    }

    pub fn get_grid(&self) -> Grid {
        self.mac_grid.grid
    }
}