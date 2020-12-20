pub mod particle;
pub mod mac_grid;
pub mod forces;
pub(crate) mod conjugate_gradient;
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
                grid_shape: (2, 2, 2),
                offset: na::Vector3::zeros()
            }
        }
    }
}

pub struct FluidSim<T> where T: Force {
    dt: f32,
    mac_grid: MACGrid,
    velocity: na::DVector<f32>,
    pressure: na::DVector<f32>,
    pub particles: Particles,
    external_force: T
}

/// Fluid simulator
impl<T> FluidSim <T> where T: Force {

    /// #Arguments:
    /// * `config` - The configuration for the simulation.
    /// * `external_force` - The external forces acting on the fluid.
    /// # Returns:
    /// * A newly constructed fluid sim
    pub fn new(config: SimConfig, external_force: T) -> Self {
        let mac_grid = MACGrid::new(&config.grid, config.density);
        let velocity = na::DVector::zeros(mac_grid.get_velocity_size());
        let pressure = na::DVector::zeros(mac_grid.get_pressure_size());
        Self {
            dt: config.delta_time,
            mac_grid,
            particles: Particles::filled(&config.grid),
            external_force,
            velocity,
            pressure
        }
    }

    /// Results:
    /// * updates the state of the simulation by one time step of length dt.
    pub fn step(&mut self) {
        let div = self.mac_grid.div_velocity_operator();
        let div_grad = self.mac_grid.div_grad_pressure_operator();
        self.particles.advect(self.dt);
        self.particles.apply_force(&self.external_force, self.dt);
        self.mac_grid.set_from_particles(&mut self.velocity, &self.particles);
        self.particles.update_using_pic(&self.mac_grid, &self.velocity);
    }

    pub fn get_grid(&self) -> Grid {
        self.mac_grid.base_grid
    }
}