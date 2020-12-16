pub mod particle;
pub mod mac_grid;
use nalgebra as na;

pub struct FluidSim<T> {
    max_grid: mac_grid::MACGrid,
    particles: Vec<particle::Particle<T>>
}

/// Fluid simulator
impl<T> FluidSim<T> {
    /// Arguments:
    /// * dims - the dimension of the sim volume (depth, height, width). Note
    /// that all sim volumes start from 0, 0, 0.
    /// * grid_shape - This is the number of cells along each axis.
    /// 
    pub fn new(dims: na::Vector3<f32>, grid_shape: (usize, usize, usize)) -> Self {
        unimplemented!()
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
    /// Restults:
    /// * updates the state of the simulation by one time step of length dt.
    pub fn step(&mut self, dt: f32) {
        unimplemented!();
    }
}