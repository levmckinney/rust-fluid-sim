pub mod particle;
pub mod mac_grid;
pub mod forces;
pub(crate) mod conjugate_gradient;
pub mod utils;

use particle::{Particle, Particles};
use conjugate_gradient::conjugate_gradient;
use forces::Force;
use mac_grid::MACGrid;
use utils::Grid;
use nalgebra as na;

/// density - This is the number of cells along each axis. 
/// dims - The dimension of the sim volume (depth, height, width).
/// Note that all sim volumes start from 0, 0, 0.
#[derive(Clone, Copy)]
pub struct SimConfig {
    delta_time: f64,
    density: f64,
    grid: Grid
}

impl Default for SimConfig {
    fn default() -> Self {
        SimConfig {
            delta_time: 1.0/60.0,
            density: 1.225, // density of air in Kg per meter cubed
            grid: Grid {
                dims: na::Vector3::new(1.0, 1.0, 1.0),
                grid_shape: (11, 11, 11),
                offset: na::Vector3::zeros()
            }
        }
    }
}

pub struct FluidSim<T> where T: Force {
    dt: f64,
    density: f64,
    pub mac_grid: MACGrid,
    pub velocity: na::DVector<f64>,
    pub pressure: na::DVector<f64>,
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
            density: config.density,
            mac_grid,
            particles: Particles::filled(&config.grid),
            external_force,
            velocity,
            pressure,
        }
    }

    /// Results:
    /// * updates the state of the simulation by one time step of length dt.
    pub fn step(&mut self) {
        let old_velocity_vec = self.velocity.clone();
        self.particles.advect(self.dt);
        self.particles.apply_force(&self.external_force, self.dt);
        self.mac_grid.set_from_particles(&mut self.velocity, &self.particles);
        println!("velocity {}", self.velocity.norm());
        let d = (self.density/self.dt)*self.mac_grid.div_velocity_operator(&self.velocity);
        println!("div velocity {}", d.norm());
        let new_pressure = conjugate_gradient(
            &self.mac_grid,
            &self.pressure,
            &d,
            1.0e-12,
            10000);
        
        self.pressure = new_pressure;
        println!("max velocity {}", self.velocity.abs().max());
        self.velocity -= (self.dt/self.density)*self.mac_grid.grad_pressure(&self.pressure);
        self.particles.update_using_flip(&self.mac_grid, &self.velocity, &old_velocity_vec);
    }

    pub fn get_grid(&self) -> Grid {
        self.mac_grid.base_grid
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use forces::ConstantForce;
    #[test]
    fn test_sum_duplicate_triplets() {
        let mut fluid_sim = FluidSim::<ConstantForce>::new(SimConfig::default(), ConstantForce{ force: na::Vector3::new(0.0, -0.00098, 0.0)});
        fluid_sim.step();
        fluid_sim.step();
        fluid_sim.step();
        fluid_sim.step();
    }
}