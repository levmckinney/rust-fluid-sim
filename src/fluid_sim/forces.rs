use super::particle::Particle;
use nalgebra as na;

pub trait Force {
    fn get_force(&self, particle: &Particle) -> na::Vector3<f64>;
}

pub struct ConstantForce {
    pub force: na::Vector3<f64>
}

impl Force for ConstantForce {
    fn get_force(&self, particle: &Particle) -> na::Vector3<f64> {
        self.force
    }
}


pub struct SphereForce {
    pub center: na::Vector3<f64>,
    pub radius: f64,
    pub inside_force: na::Vector3<f64>,
    pub outside_force: na::Vector3<f64>
}


impl Force for SphereForce {
    fn get_force(&self, particle: &Particle) -> na::Vector3<f64> {
        if (particle.position - self.center).norm().abs() < self.radius {
            self.inside_force
        } else {
            self.outside_force
        }
    }
}