use super::particle::Particle;
use nalgebra as na;

pub trait Force {
    fn get_force(&self, particle: &Particle) -> na::Vector3<f32>;
}

pub struct ConstantForce {
    pub force: na::Vector3<f32>
}

impl Force for ConstantForce {
    fn get_force(&self, particle: &Particle) -> na::Vector3<f32> {
        self.force
    }
}