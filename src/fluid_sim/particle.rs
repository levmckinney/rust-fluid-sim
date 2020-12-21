use super::mac_grid::{MACGrid, Component};
use super::forces::Force;
use super::utils::{Grid};

use rand::{thread_rng, Rng};
use rand::distributions::Uniform;
use rand::distributions::Distribution;
use nalgebra as na;

#[derive(Clone, Copy)]
pub struct Particle {
    pub position: na::Vector3<f64>,
    pub velocity: na::Vector3<f64>
}

/// Store particles on a sub grid which divides
/// each cell of the simulation grid into 8 sub cells
pub struct Particles {
    pub sub_grid: Grid,
    particles: Vec<Vec<Vec<Vec<Particle>>>>
}


impl Particles {
    // Public
    /// Creates empty grid of particles
    pub fn new(grid: &Grid) -> Self {
        let sub_grid = Grid {
            grid_shape: (2*grid.grid_shape.0, 2*grid.grid_shape.1, 2*grid.grid_shape.1),
            ..*grid
        };
        let mut particles = Vec::new();
        particles.reserve(sub_grid.grid_shape.0);
        for _ in 0..sub_grid.grid_shape.0 {            
            let mut layer = Vec::new();
            particles.reserve(sub_grid.grid_shape.0);
            for _ in 0..sub_grid.grid_shape.1 {
                let mut col = Vec::new();
                for _ in 0..sub_grid.grid_shape.2 {
                    let cell = Vec::new();
                    col.push(cell);
                }
                layer.push(col);
            }
            particles.push(layer);
        }

        Particles {
            sub_grid: sub_grid, 
            particles: particles
        }
    }

    /// Create particles within grid cells in the same pattern as described within 
    /// the flip paper i.e. every cell has 8 particles which are randomly positioned within
    /// their respective sub cells.
    pub fn filled(grid: &Grid) -> Self {
        let mut particles = Particles::new(&grid);
        let cell_size = grid.cell_size();
        let mut rng = thread_rng();
        let uniform = Uniform::new(0.0, 1.0);
        for (i, j, k) in particles.sub_grid.cells() {
            let sub_cell_pos = particles.sub_grid.get_position(i, j, k);
            let jitter = na::Vector3::<f64>::new(
                uniform.sample(&mut rng)*cell_size[0],
                uniform.sample(&mut rng)*cell_size[1],
                uniform.sample(&mut rng)*cell_size[2]
            );
            particles.particles[i][j][k].push(Particle {
                position: sub_cell_pos + jitter,
                velocity: na::Vector3::<f64>::new(0.0, 0.0, 0.0),
            })
        }
        particles
    }
    
    /// Update the particles positions using forward euler.
    pub fn advect(&mut self, dt: f64) {
        for (i, j, k) in self.sub_grid.cells() {
            for particle in &mut self.particles[i][j][k] {
                particle.position += particle.velocity*dt;
            }
        }
        self.fix_particles();
    }

    /// # Arguments
    /// * mac_grid - a mac grid representing static grid.
    /// * q - the generalized velocities of the grid
    /// # Results
    /// * Updates the particles velocities using the pic method.
    pub fn update_using_pic(&mut self, mac_grid: &MACGrid, velocity_vec: &na::DVector<f64>) {
        for (i, j, k) in self.sub_grid.cells() {
            for comp in Component::iterator() {
                let comp_grid = mac_grid.get_velocity_grid(comp);
                let enclosing_grid_cell = match comp {
                    Component::U => ((i + 1)/2, j/2, k/2),
                    Component::V => (i/2, (j + 1)/2, k/2),
                    Component::W => (i/2, j/2, (k + 1)/2)
                };
                
                let corners = iproduct!(0..2, 0..2, 0..2).map(|c|
                    (enclosing_grid_cell.0 + c.0, enclosing_grid_cell.1 + c.1, enclosing_grid_cell.2 + c.2)
                ).collect::<Vec<(usize, usize, usize)>>();
                // Now we linearly interpolate from the corners to find the velocity at the particles position
                // within the grid cell.
                for particle in &mut self.particles[i][j][k] {
                    let mut velocity = 0.0;
                    let mut sum = 0.0;
                    for corner in &corners {
                        let weight = comp_grid.bilinear_weight(corner, &particle.position);
                        let vel_component = velocity_vec[mac_grid.get_velocity_ind(comp, corner.0, corner.1, corner.2)];
                        sum += weight;
                        velocity += weight*vel_component
                    }
                    velocity /= sum;
                    match comp {
                        Component::U => {
                            particle.velocity[0] = velocity;
                        }
                        Component::V => {
                            particle.velocity[1] = velocity;
                        }
                        Component::W => {
                            particle.velocity[2] = velocity;
                        }
                    }
                }
            }
        }
    }

    /// # Arguments
    /// * `force` - force that can be sampled at each particle
    /// Apply external force to all particles using forward euler.
    pub fn apply_force<T>(&mut self, force: &T, dt: f64) where T: Force{
        for (i, j, k) in self.sub_grid.cells() {
            for particle in &mut self.particles[i][j][k] {
                particle.velocity += dt*force.get_force(particle);
            }
        }
    }

    /// # Arguments:
    /// * `min` - minimum extent of the sub grid we want to sample
    /// * `max` - maximum extent of the sub grid we want to sample
    /// # Returns:
    /// * A Vector contain all particles within the grid bounded by max and min inclusive
    pub fn get_particles_within(&self, min: (usize, usize, usize), max: (usize, usize, usize)) -> Vec::<&Particle> {
        let min = self.sub_grid.clamp_inds(min.0, min.1, min.2);
        let max = self.sub_grid.clamp_inds(max.0, max.1, max.2);

        let mut particles_within = Vec::new();

        for (i, j, k) in iproduct!(min.0..(max.0 + 1), min.1..(max.1 + 1), min.2..(max.2 + 1)) {
            // TODO remove extra assert
            for particle in &self.particles[i][j][k] {
                debug_assert_eq!(self.sub_grid.get_containing_cell_ind(&particle.position), (i, j, k), "Particle out of proper sub grid cell!")
            }
            particles_within.extend(&self.particles[i][j][k]);
        }
        particles_within
    }


    /// # Returns
    /// * A vector containing a copy of every particle in the particle grid
    pub fn get_particles(&self) -> Vec::<Particle> {
        let mut particle_vec = Vec::new();
        for (i,j,k) in self.sub_grid.cells() {
            particle_vec.extend((*self.particles[i][j][k]).into_iter().cloned());
        }
        particle_vec
    }

    // Private
    // Ensure that all particles are within there respective grid cells
    fn fix_particles(&mut self) {
        let mut out_of_cells = Vec::<Particle>::new();
        for (i, j, k) in self.sub_grid.cells() {
            let corner1 = self.sub_grid.get_position(i, j, k);
            let corner2 = self.sub_grid.get_position(i+1, j+1, k+1);
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
        // Place particles that lay outside of their cells into their place.
        for particle in out_of_cells {
            let (i, j, k) = self.sub_grid.get_containing_cell_ind(&particle.position);
            self.particles[i][j][k].push(
                Particle {position: self.sub_grid.clamp(&particle.position), velocity: particle.velocity}
            );
        }
    }
}

#[cfg(test)]
mod tests {

}
