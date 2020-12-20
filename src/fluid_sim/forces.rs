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


pub struct SphereForce {
    pub center: na::Vector3<f32>,
    pub radius: f32,
    pub force: na::Vector3<f32>,
}


impl Force for SphereForce {
    fn get_force(&self, particle: &Particle) -> na::Vector3<f32> {
        if (particle.position - self.center).norm().abs() < self.radius {
            self.force
        } else {
            na::Vector3::zeros()
        }
    }
}