use super::mac_grid::{MACGrid};
use super::forces::Force;
use super::utils::{Grid, AABB, clamp};

use std::iter;
use rand::{thread_rng, Rng};
use rand::distributions::Uniform;
use rand::distributions::Distribution;
use nalgebra as na;

#[derive(Clone, Copy)]
pub struct Particle {
    pub position: na::Vector3<f32>,
    pub velocity: na::Vector3<f32>
}

pub struct Particles {
    pub grid: Grid,
    pub particles: Vec<Vec<Vec<Vec<Particle>>>>
}


impl Particles {
    pub fn new(grid: &Grid) -> Self {
        let mut particles = Vec::new();
        particles.reserve(grid.grid_shape.0);
        for _ in 0..grid.grid_shape.0 {            
            let mut layer = Vec::new();
            particles.reserve(grid.grid_shape.0);
            for _ in 0..grid.grid_shape.1 {
                let mut col = Vec::new();
                for _ in 0..grid.grid_shape.2 {
                    let cell = Vec::new();
                    col.push(cell);
                }
                layer.push(col);
            }
            particles.push(layer);
        }

        Particles {
            grid: *grid, 
            particles: particles
        }
    }

    /// Create particles within grid cells in the same pattern as described within 
    /// the flip paper i.e. every cell has 8 particles which are randomly positioned within
    /// their respective sub cells.
    pub fn filled(grid: &Grid) -> Self {
        let Grid {dims, grid_shape} = grid;
        let cell_size = na::Vector3::<f32>::new(
            dims[0]/(grid_shape.0 as f32), 
            dims[1]/(grid_shape.1 as f32), 
            dims[2]/(grid_shape.2 as f32));
        let sub_cell_size = cell_size/2.0;
        let mut particles = Particles::new(&grid);
        let mut rng = thread_rng();
        let uniform = Uniform::new(0.0, 1.0);
        for (i, j, k) in grid.cells() {
            let cell_pos = na::Vector3::<f32>::new(
                cell_size[0]*(i as f32),
                cell_size[1]*(j as f32),
                cell_size[2]*(k as f32));
            for (i_sub, j_sub, k_sub) in iproduct!(0..2, 0..2, 0..2) {
                let sub_cell_pos = na::Vector3::<f32>::new(
                    sub_cell_size[0]*(i_sub as f32),
                    sub_cell_size[1]*(j_sub as f32),
                    sub_cell_size[2]*(k_sub as f32));
                let jitter = na::Vector3::<f32>::new(
                    uniform.sample(&mut rng)*sub_cell_size[0],
                    uniform.sample(&mut rng)*sub_cell_size[1],
                    uniform.sample(&mut rng)*sub_cell_size[2]
                );
                particles.particles[i][j][k].push(Particle {
                    position: cell_pos + sub_cell_pos + jitter,
                    velocity: na::Vector3::<f32>::new(0.0, 0.0, 0.0),
                })
            }
        }
        particles
    }

    // Ensure that all particles are within there respective grid cells
    fn fix_particles(&mut self) {
        let mut out_of_cells = Vec::<Particle>::new();
        let cell_size = self.grid.cell_size();
        for (i, j, k) in self.grid.cells() {
            let corner1 = na::Vector3::<f32>::new(
                cell_size[0]*(i as f32),
                cell_size[1]*(j as f32),
                cell_size[2]*(k as f32));
            let corner2 = na::Vector3::<f32>::new(
                cell_size[0]*((i + 1) as f32),
                cell_size[1]*((j + 1) as f32),
                cell_size[2]*((k + 1) as f32));
            // Split the particles in this cell into particles that are supposed to be in this cell and those that are not.
            let mut in_cell= Vec::new();
            for particle in &self.particles[i][j][k] {
                if corner1 < particle.position && corner2 > particle.position {
                    in_cell.push(*particle)
                } else {
                    out_of_cells.push(*particle);
                }    
            }
            self.particles[i][j][k] = in_cell;
        }
        // Place particles that lay outside of there cells into their place.
        let zero =  na::Vector3::<f32>::new(0.0, 0.0, 0.0);
        for particle in out_of_cells {
            let clamped_position = clamp(particle.position, zero, self.grid.dims);
            let (i, j, k) = (
                (clamped_position[0] / cell_size[0]) as usize,
                (clamped_position[1] / cell_size[1]) as usize,
                (clamped_position[2] / cell_size[2]) as usize);
            assert!(self.grid.within(i, j, k));
            self.particles[i][j][k].push(
                Particle {position: clamped_position, velocity: particle.velocity}
            );
        }
    }

    /// Update the particles positions using advection
    pub fn advect(&mut self, dt: f32) {
        for (i, j, k) in self.grid.cells() {
            for particle in &mut self.particles[i][j][k] {
                particle.position += particle.velocity;
            }
        }

        self.fix_particles();
    }

    /// Update the particles velocities using flip.
    pub fn update_using_flip(&mut self, mac_grid: &MACGrid) {
        unimplemented!();
    }

    /// Apply external force to all particles using forward euler.
    pub fn apply_force<T>(&mut self, force: &T, dt: f32) where T: Force{
        for (i, j, k) in self.grid.cells() {
            for particle in &mut self.particles[i][j][k] {
                particle.velocity += dt*force.get_force(particle);
            }
        }
    }

    pub fn get_particles_within(aabb: AABB) -> Vec<Particle> {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_vec_comp() {
        assert_eq!(na::Vector3::new(0.0, 0.0, 0.0) < na::Vector3::new(1.0, 1.0, 1.0), true)
    }
}
