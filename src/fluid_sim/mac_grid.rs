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
    pub fn get_velocity_size(&self) -> usize {
        self.u_grid.num_cells() + self.v_grid.num_cells() + self.w_grid.num_cells()
    }

    /// Get the size of the pressure vector required by this mac_grid.
    pub fn get_pressure_size(&self) -> usize {
        self.p_grid.num_cells()
    }

    // Public interface
    // Getters and setters
    pub fn get_velocity_grid(&self, c: &Component) -> &Grid {
        match c {
            Component::U => &self.u_grid,
            Component::V => &self.v_grid,
            Component::W => &self.w_grid
        }
    }
    /// Get the index for component_ijk in the generalized velocity vector
    pub fn get_velocity_ind(&self, c: &Component, i: usize, j: usize, k: usize) -> usize {
        match c {
            Component::U => self.u_grid.flat_ind(i, j, k),
            Component::V => self.v_grid.flat_ind(i, j, k) + self.u_grid.num_cells(),
            Component::W => self.w_grid.flat_ind(i, j, k) + self.u_grid.num_cells() + self.v_grid.num_cells()
        }
    }

    /// Get the index for pressure at p_ijk in the pressure vector
    pub fn get_pressure_ind(&self, i: usize, j: usize, k: usize) -> usize {
        self.p_grid.flat_ind(i, j, k)
    }

    pub fn set_from_particles(&self, q: &mut na::DVector<f32>, particles: &Particles){
        for comp in Component::iterator() {
            let comp_grid = self.get_velocity_grid(comp);
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
                q[self.get_velocity_ind(comp, i, j, k)] = velocity;
            }
        }
    }

    /// Assembles the global divergence operator matrix.
    pub fn div_velocity_operator(&self) -> na::CsMatrix<f32> {
        let div_op_local = self.div_operator_local();
        let mut triplets = Vec::new();
        for (row, (_i, _j, _k)) in self.p_grid.cells().enumerate() {
            if self.is_boundary_pressure(_i, _j, _k) {
                continue;
            }
            // Move from working with indices over the pressure grid to
            // indices over the base grid.
            let (i, j, k) = (_i-1, _j-1, _k-1);
            let boundary_filter = self.get_boundary_filter(i, j, k);
            let velocity_grid_inds = self.get_local_velocity_inds(i, j, k);
            // Convert the indices over velocity grids to indices
            // into the generalized velocity vector.
            let velocity_vector_inds = velocity_grid_inds.iter().map(|(comp, i, j, k)| 
                self.get_velocity_ind(comp, *i, *j, *k));

            let operator = div_op_local*boundary_filter;

            triplets.extend(
                velocity_vector_inds.zip(operator.iter()).map(|(vel_ind, value)|
                    ((row, vel_ind), *value)));
        }
        sum_from_triplets(self.get_pressure_size(), self.get_velocity_size(), triplets)
    }

    /// Returns global divergence of grad of pressure operator.
    pub fn div_grad_pressure_operator(&self) -> na::CsMatrix<f32>{
        let div_op_local = self.div_operator_local();
        let grad_op_local = self.grad_operator_local();
        let mut triplets = Vec::new();
        for (row, (i, j, k)) in self.p_grid.cells().enumerate() {
            if self.is_boundary_pressure(i, j, k) {
                continue;
            }
            let pressure_grid_inds = self.get_local_pressure_inds(i, j, k);
            // Subtract one here because get_boundary_filter takes inds over 
            // the base grid.
            let boundary_filter = self.get_boundary_filter(i - 1, j - 1, k - 1);
            // Convert indices over the pressure grid to indices over the
            // the pressure vector
            let pressure_vector_inds = pressure_grid_inds.iter().map(|(i, j, k)| 
                self.get_pressure_ind(*i, *j, *k));

            let operator = div_op_local*boundary_filter*grad_op_local;

            triplets.extend(
                pressure_vector_inds.zip(operator.into_iter()).map(|(pressure_ind, value)|
                    ((row, pressure_ind), *value)));
        }
        sum_from_triplets(self.get_pressure_size(), self.get_pressure_size(), triplets)
    }
    
    //Private methods
    /// Returns weather velocity component c at position i, j, k is on or outside the boundary.
    fn is_boundary_velocity(&self, c: &Component, i: usize, j: usize, k: usize) -> bool {
        let comp_grid = self.get_velocity_grid(c);
        (i == 0) 
        || (j == 0) 
        || (k == 0) 
        || (i >= comp_grid.grid_shape.0 - 1) 
        || (j >= comp_grid.grid_shape.1 - 1) 
        || (k >= comp_grid.grid_shape.2 - 1)
    }

    // Useful latter for implementing constant flow boundary condition.
    fn boundary_velocity(&self, c: &Component, i: usize, j:usize, k: usize) -> f32 {
        0.0
    }

    /// Returns true pressure at point i, j, k is outside or on the boundary.
    fn is_boundary_pressure(&self, i: usize, j:usize, k: usize) -> bool {
        (i == 0) 
        || (j == 0) 
        || (k == 0) 
        || (i >= self.p_grid.grid_shape.0 - 1)
        || (j >= self.p_grid.grid_shape.1 - 1)
        || (k >= self.p_grid.grid_shape.2 - 1)
    }
    
    /// Take a base cell and return the local velocity indices. 
    fn get_local_velocity_inds(&self, i: usize, j: usize, k: usize) -> [(Component, usize, usize, usize); 6] {
        debug_assert!(self.base_grid.index_within(i, j, k));
        [(Component::U, i + 1, j, k), 
        (Component::U, i + 2, j, k), 
        (Component::V, i, j + 1, k), 
        (Component::V, i, j + 2, k), 
        (Component::W, i, j, k + 1),
        (Component::W, i, j, k + 2)]
    }

    /// Take an index into pressures and return the indices of the surounding pressures.
    fn get_local_pressure_inds(&self, i: usize, j: usize, k: usize) -> [(usize, usize, usize); 7] {
        debug_assert!(!self.is_boundary_pressure(i, j, k));
        // Index in the pressure grid
        [(i, j, k), 
         (i - 1, j, k),
         (i + 1, j, k),
         (i, j - 1, k),
         (i, j + 1, k),
         (i, j, k - 1),
         (i, j, k + 1)]
    }

    /// Returns the boundary filter for the cell at position 
    fn get_boundary_filter(&self, i: usize, j: usize, k: usize) -> na::Matrix6<f32> {
        let not_on_boundary = na::Vector6::<f32>::from_iterator(
            self.get_local_velocity_inds(i, j, k).iter().map(|(comp, _i, _j, _k)|
                if self.is_boundary_velocity(comp, *_i, *_j, *_k) {1.0} else {0.0}));
        na::Matrix6::from_diagonal(&not_on_boundary)
    }

    /// Returns local div operator
    fn div_operator_local(&self) -> na::Matrix1x6<f32> {
        let delta = self.base_grid.cell_size();
        na::Matrix1x6::new(
            -1.0/delta[0], 1.0/delta[0], -1.0/delta[1], 1.0/delta[1], -1.0/delta[2], 1.0/delta[2] 
        )
    }

    // return the local divergence of gradient operator
    fn grad_operator_local(&self) -> na::MatrixMN<f32, U6, U7> {
        let delta = self.base_grid.cell_size();
        let mut grad_op = na::MatrixMN::<f32, U6, U7>::zeros();
        grad_op[(0, 0)] = 1.0/delta[0]; grad_op[(0, 1)] = -1.0/delta[0];
        grad_op[(1, 0)] = -1.0/delta[0]; grad_op[(1, 2)] = 1.0/delta[0];
        grad_op[(2, 0)] = 1.0/delta[1]; grad_op[(2, 3)] = -1.0/delta[1];
        grad_op[(3, 0)] = -1.0/delta[1]; grad_op[(3, 4)] = 1.0/delta[1];
        grad_op[(4, 0)] = 1.0/delta[2]; grad_op[(4, 5)] = -1.0/delta[2];
        grad_op[(5, 0)] = -1.0/delta[2]; grad_op[(5, 6)] = 1.0/delta[2];
        grad_op
    }
}

