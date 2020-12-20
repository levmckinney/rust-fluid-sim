use super::particle::{Particle, Particles};
use super::utils::{sum_from_triplets, Triplets, Grid};
use nalgebra::base::dimension::{U6, U7};
use std::slice::Iter;
use nalgebra as na;

#[derive(Debug)]
pub enum Component {
    U,
    V,
    W
}

impl Component {
    pub fn iterator() -> Iter<'static, Component> {
        static COMPONENT: [Component; 3] = [Component::U, Component::V, Component::W];
        COMPONENT.iter()
    }
}

pub struct MACGrid {
    pub base_grid: Grid,
    u_grid: Grid,
    v_grid: Grid,
    w_grid: Grid,
    p_grid: Grid,
    pub density: f32,
}

impl MACGrid {
    pub fn new(grid: &Grid, density: f32) -> Self {
        let Grid {dims, grid_shape, offset} = grid;
        let cell_size = grid.cell_size();
        let u_grid = Grid::new(
            offset - na::Vector3::new(cell_size[0]/2.0, 0.0, 0.0),
            cell_size,
            (grid_shape.0 + 2, grid_shape.1 + 1, grid_shape.2 + 1)
        );
        let v_grid = Grid::new(
            offset - na::Vector3::new(0.0, cell_size[1]/2.0, 0.0),
            cell_size,
            (grid_shape.0 + 1, grid_shape.1 + 2, grid_shape.2 + 1)
        );
        let w_grid = Grid::new(
            offset - na::Vector3::new(0.0, 0.0, cell_size[2]/2.0),
            cell_size,
            (grid_shape.0 + 1, grid_shape.1 + 1, grid_shape.2 + 2)
        );
        let p_grid = Grid::new(
            offset - cell_size*0.5,
            cell_size,
            (grid_shape.0 + 2, grid_shape.1 + 2, grid_shape.2 + 2)
        );

        MACGrid {
            base_grid: *grid,
            u_grid,
            v_grid,
            w_grid, 
            p_grid,
            density
        }
    }

    /// Get the size of the generalized velocity vector 
    /// required by this mac_grid
    pub fn get_q_size(&self) -> usize {
        self.u_grid.num_cells() + self.v_grid.num_cells() + self.w_grid.num_cells()
    }

    /// Get the size of the pressure vector required by this mac_grid.
    pub fn get_p_size(&self) -> usize {
        self.w_grid.num_cells()
    }

    // Public interface
    // Getters and setters
    pub fn get_vel_grid(&self, c: &Component) -> &Grid {
        match c {
            Component::U => &self.u_grid,
            Component::V => &self.v_grid,
            Component::W => &self.w_grid
        }
    }
    /// Get the index for component_ijk in the generalized coordinates q
    pub fn get_vel_ind(&self, c: &Component, i: usize, j: usize, k: usize) -> usize {
        match c {
            Component::U => self.u_grid.flat_ind(i, j, k),
            Component::V => self.v_grid.flat_ind(i, j, k) + self.u_grid.num_cells(),
            Component::W => self.w_grid.flat_ind(i, j, k) + self.u_grid.num_cells() + self.v_grid.num_cells()
        }
    }

    pub fn set_from_particles(&self, q: &mut na::DVector<f32>, particles: &Particles){
        for comp in Component::iterator() {
            let comp_grid = self.get_vel_grid(comp);
            for (i, j, k) in comp_grid.cells() {
                let mut velocity = 0.0;
                if self.is_boundary_velocity(comp, i, j, k) {
                    velocity = self.boundary_velocity(comp, i, j, k);
                } else {
                    // The index of the corresponding cell in the
                    // particle subdivided grid.
                    let sub_grid_ind = match comp {
                        Component::U => (2*i - 1, 2*j, 2*k),
                        Component::V => (2*i, 2*j - 1, 2*k),
                        Component::W => (2*i, 2*j, 2*k - 1)
                    }; 
                    // We get all particles within a 2 by 2 grid cell box surrounding 
                    // the position of the velocity component
                    let particles_nearby = particles.get_particles_within(
                        // We need to use saturating subtract to avoid underflow.
                        (sub_grid_ind.0.saturating_sub(2), sub_grid_ind.1.saturating_sub(2), sub_grid_ind.2.saturating_sub(2)),
                        // Going out of bounds on the maximum is well defined by get particles is.
                        (sub_grid_ind.0 + 1, sub_grid_ind.1 + 1, sub_grid_ind.2 + 1)
                    );
                    // Get weighted average of nearby particles velocity components.
                    let mut sum = 0.0;
                    for particle in particles_nearby {
                        let particle_velocity_comp = match comp {
                            Component::U => particle.velocity[0],
                            Component::V => particle.velocity[1],
                            Component::W => particle.velocity[2],
                        };
                        let weight = comp_grid.bilinear_weight(&(i, j, k), &particle.position);
                        sum += weight;
                        velocity += weight*particle_velocity_comp;
                    }
                    velocity /= sum;
                }
                // Set the proper index in velocity to newly interpolated value.
                q[self.get_vel_ind(comp, i, j, k)] = velocity;
            }
        }
    }

    pub fn pressure_project(&self) {
        unimplemented!();
    }

    //Private methods
    /// Returns weather velocity component c at position i, j, k is on or outside the boundary.
    fn is_boundary_velocity(&self, c: &Component, i: usize, j: usize, k: usize) -> bool {
        match c {
            Component::U => (i == 0) 
                || (j == 0) 
                || (k == 0) 
                || (i >= self.base_grid.grid_shape.0 + 1) 
                || (j >= self.base_grid.grid_shape.1) 
                || (k >= self.base_grid.grid_shape.2),
            Component::V => (i == 0) 
                || (j == 0) 
                || (k == 0) 
                || (i >= self.base_grid.grid_shape.0) 
                || (j >= self.base_grid.grid_shape.1 + 1) 
                || (k >= self.base_grid.grid_shape.2),
            Component::W => (i == 0) 
                || (j == 0) 
                || (k == 0) 
                || (i >= self.base_grid.grid_shape.0) 
                || (j >= self.base_grid.grid_shape.1) 
                || (k >= self.base_grid.grid_shape.2 + 1)
        }
        
    }

    /// Returns true pressure at point i, j, k is outside or on the boundary.
    fn is_boundary_pressure(&self, i: usize, j:usize, k: usize) -> bool {
        unimplemented!();
    }

    // Useful latter for implementing constant flow boundary condition.
    fn boundary_velocity(&self, c: &Component, i: usize, j:usize, k: usize) -> f32 {
        0.0
    }


    /// Returns the boundary filter for the cell at position 
    fn get_boundary_filter(&self, i:usize, j:usize, k:usize) -> na::Matrix6<f32> {
        unimplemented!();
    }

    /// Returns local div operator
    fn div_operator_local(&self, deltas: na::Vector3<f32>) -> na::Matrix6x1<f32> {
        unimplemented!();
    }

    // what is the shape of the input operator vector?
    fn div_of_grad_operator_local(&self, deltas: na::Vector3<f32>) -> na::MatrixMN<f32, U6, U7> {
        unimplemented!();
    }

    /// Returns global divergence operator matrix.
    fn div_of_velocity_operator_global(&self, deltas: na::Vector3<f32>) -> na::CsMatrix<f32> {
        unimplemented!();
    }

    /// Returns global divergence of grad of pressure matrix.
    fn div_of_grad_of_pressure_operator_global(&self, deltas: na::Vector3<f32>) -> na::CsMatrix<f32>{
        unimplemented!();
    }
}

