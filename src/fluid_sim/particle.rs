use nalgebra as na;
use super::mac_grid::{MACGrid};

pub struct Particle<T> {
    pub params: T,
    pub position: na::Vector3<f32>,
    pub velocity: na::Vector3<f32>
}

pub struct Particles<T> {
    pub particles: Vec<Particle<T>>
}

impl<T> Particles<T> {
    /// Update the particles positions using advection
    pub fn advect(&mut self) {
        unimplemented!();
    }

    /// Update the particles velocities using flip.
    pub fn update_using_flip(&mut self, mac_grid: &MACGrid) {
        unimplemented!();
    }
}
