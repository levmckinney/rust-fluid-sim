use super::particle::{Particles};
use super::utils::{Grid};
use std::slice::Iter;

use nalgebra::base::dimension::{U6, U7};
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

/// A struct representing a mac grid
pub struct MACGrid {
    /// The space of cells to simulate.
    pub base_grid: Grid, 
    /// Staggered grid velocities in the X direction.
    pub u_grid: Grid, 
    /// Staggered grid velocities in the Y direction.
    pub v_grid: Grid, 
    /// Staggered grid velocities in the Z direction.
    pub w_grid: Grid, 
    /// Staggered grid holding pressures at the center of each cell.
    pub p_grid: Grid, 
}

impl MACGrid {
    // Public interface
    pub fn new(grid: &Grid) -> Self {
        let Grid {dims:_, grid_shape, offset} = grid;
        let cell_size = grid.cell_size();
        // The staggered grids
        let u_grid = Grid::new(
            offset - na::Vector3::new(0.0, cell_size[0]/2.0, cell_size[0]/2.0),
            cell_size,
            (grid_shape.0 + 1, grid_shape.1 + 2, grid_shape.2 + 2)
        );
        let v_grid = Grid::new(
            offset - na::Vector3::new(cell_size[0]/2.0, 0.0, cell_size[0]/2.0),
            cell_size,
            (grid_shape.0 + 2, grid_shape.1 + 1, grid_shape.2 + 2)
        );
        let w_grid = Grid::new(
            offset - na::Vector3::new(cell_size[0]/2.0, cell_size[0]/2.0, 0.0),
            cell_size,
            (grid_shape.0 + 2, grid_shape.1 + 2, grid_shape.2 + 1)
        );
        // The pressure grid is 2 larger on in all directions to
        // account for the pressures within the walls.
        // Note that the pressure grid is offset buy one cell. So point (0,0,0) in
        // the base grid is (1, 1, 1) in the pressure grid.
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
        }
    }

    /// Get the size of the generalized velocity vector required by this mac_grid
    pub fn get_velocity_size(&self) -> usize {
        self.u_grid.num_cells() + self.v_grid.num_cells() + self.w_grid.num_cells()
    }

    /// Get the size of the pressure vector required by this mac_grid.
    pub fn get_pressure_size(&self) -> usize {
        self.p_grid.num_cells()
    }

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

    /// Set each component the velocity grid based on the 
    pub fn set_from_particles(&self, velocity_vec: &mut na::DVector<f64>, particles: &Particles){
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
                        Component::U => (2*i, 2*j - 1, 2*k - 1),
                        Component::V => (2*i - 1, 2*j, 2*k - 1),
                        Component::W => (2*i - 1, 2*j - 1, 2*k)
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
                    let mut sum = 0.001;
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
                velocity_vec[self.get_velocity_ind(comp, i, j, k)] = velocity;
            }
        }
    }

    /// # Arguments
    /// * `velocity_vec` - a self.get_velocity_size() vector representing the velocity at every point within the grid and
    ///     on the boundary
    /// # Returns
    ///     the velocity divergence of the velocity
    pub fn div_velocity_operator(&self, velocity_vec: &na::DVector<f64>) -> na::DVector<f64> {
        let div_op_local = self.div_operator_local();
        na::DVector::from_iterator(self.get_pressure_size(), 
            self.p_grid.cells().map(|(_i, _j, _k)| {
                if self.is_boundary_pressure(_i, _j, _k) {
                    return 0.0;
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
                let local_velocities = na::Vector6::<f64>::from_iterator(
                    velocity_vector_inds.map(|l| velocity_vec[l]));
                (div_op_local*boundary_filter*local_velocities)[0]}))
    }

    /// # Arguments
    /// * `pressure_vec` - a self.get_pressure_size() dimensional vector representing
    ///     the pressure on the grid.
    /// Returns
    ///     The divergence of the gradient of the pressure.
    pub fn div_grad_pressure(&self, pressure_vec: &na::DVector<f64>) -> na::DVector<f64> {
        let div_op_local = self.div_operator_local();
        let grad_op_local = self.grad_operator_local();
        na::DVector::from_iterator(self.p_grid.num_cells(), 
            self.p_grid.cells().map(|(i, j, k)| {
                if self.is_boundary_pressure(i, j, k) {
                    return 0.0;
                }
                let boundary_filter = self.get_boundary_filter(i - 1, j - 1, k - 1);                
                // Create the local pressure vector
                let local_pressures = self.get_local_pressures(pressure_vec, i, j, k);
                // Return the result of multiplying this row in the sparse matrix
                (div_op_local*boundary_filter*grad_op_local*local_pressures)[0]
            }))
    }
    /// # Arguments
    /// * `pressure_vec` - a self.get_pressure_size() dimensional vector representing
    ///     the pressure on the grid.
    /// # Returns
    ///     A self.get_velocity_size() dimensional vector containing the gradient of the pressure.

    pub fn grad_pressure(&self, pressure_vec: &na::DVector<f64>) -> na::DVector<f64> {
        let grad_op_local = self.grad_operator_local();
        let mut grad_pressure = na::DVector::zeros(self.get_velocity_size());
        for (i, j, k) in self.base_grid.cells() {
            let boundary_filter = self.get_boundary_filter(i, j, k);
            let local_pressures = self.get_local_pressures(pressure_vec, i + 1, j + 1, k + 1) ;
            let local_grad_pressure = boundary_filter*grad_op_local*local_pressures;
            let velocity_grid_inds = self.get_local_velocity_inds(i, j, k);
            let velocity_vector_inds = velocity_grid_inds.iter().map(|(comp, i, j, k)| 
                self.get_velocity_ind(comp, *i, *j, *k));
            for (value, l) in local_grad_pressure.into_iter().zip(velocity_vector_inds) {
                // We need to divide by 2 here to avoid double counting the pressure gradient.
                // Every velocity component appearers twice in this sum because all velocities that
                // are not on the boundary are local velocities to two grid cells. 
                grad_pressure[l] += value/2.0;
            }
        }
        grad_pressure
    }

    /// Returns true pressure at point i, j, k is outside or on the boundary  of the
    /// pressure grid.
    pub fn is_boundary_pressure(&self, i: usize, j:usize, k: usize) -> bool {
        (i == 0) 
        || (j == 0) 
        || (k == 0) 
        || (i >= self.p_grid.grid_shape.0 - 1)
        || (j >= self.p_grid.grid_shape.1 - 1)
        || (k >= self.p_grid.grid_shape.2 - 1)
    }

    /// Returns weather velocity component c at position i, j, k is on or outside the boundary.
    pub fn is_boundary_velocity(&self, c: &Component, i: usize, j: usize, k: usize) -> bool {
        let comp_grid = self.get_velocity_grid(c);
        (i == 0) 
        || (j == 0) 
        || (k == 0) 
        || (i >= comp_grid.grid_shape.0 - 1) 
        || (j >= comp_grid.grid_shape.1 - 1) 
        || (k >= comp_grid.grid_shape.2 - 1)
    }
    
    //Private methods
    /// Takes indices over the pressure grid.
    /// returns local pressure vector
    fn get_local_pressures(&self, pressure_vec: &na::DVector<f64>, i: usize, j: usize, k: usize) -> na::VectorN<f64, U7> {
        let pressure_grid_inds = self.get_local_pressure_inds(i, j, k);
        let pressure_vector_inds = pressure_grid_inds.iter().map(|(i, j, k)|
            self.get_pressure_ind(*i, *j, *k));        
        // Create the local pressure vector
        na::VectorN::<f64, U7>::from_iterator(
            pressure_vector_inds.map(|l| pressure_vec[l]))
    }

    // Useful latter for implementing constant flow boundary condition.
    fn boundary_velocity(&self, _c: &Component, _i: usize, _j:usize, _k: usize) -> f64 {
        0.0
    }
    
    /// Take a base cell and return the local velocity indices. 
    fn get_local_velocity_inds(&self, i: usize, j: usize, k: usize) -> [(Component, usize, usize, usize); 6] {
        debug_assert!(self.base_grid.index_within(i, j, k));
        [(Component::U, i, j + 1, k + 1), 
        (Component::U, i + 1, j + 1, k + 1), 
        (Component::V, i + 1, j, k + 1), 
        (Component::V, i + 1, j + 1, k + 1), 
        (Component::W, i + 1, j + 1, k),
        (Component::W, i + 1, j + 1, k + 1)]
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
    fn get_boundary_filter(&self, i: usize, j: usize, k: usize) -> na::Matrix6<f64> {
        let not_on_boundary = na::Vector6::<f64>::from_iterator(
            self.get_local_velocity_inds(i, j, k).iter().map(|(comp, _i, _j, _k)|
                if self.is_boundary_velocity(comp, *_i, *_j, *_k) {0.0} else {1.0}));
        na::Matrix6::from_diagonal(&not_on_boundary)
    }

    /// Returns local div operator
    fn div_operator_local(&self) -> na::Matrix1x6<f64> {
        let delta = self.base_grid.cell_size();
        na::Matrix1x6::new(
            -1.0/delta[0], 1.0/delta[0], -1.0/delta[1], 1.0/delta[1], -1.0/delta[2], 1.0/delta[2] 
        )
    }

    // return the local divergence of gradient operator
    fn grad_operator_local(&self) -> na::MatrixMN<f64, U6, U7> {
        let delta = self.base_grid.cell_size();
        let mut grad_op = na::MatrixMN::<f64, U6, U7>::zeros();
        grad_op[(0, 0)] = 1.0/delta[0]; grad_op[(0, 1)] = -1.0/delta[0];
        grad_op[(1, 0)] = -1.0/delta[0]; grad_op[(1, 2)] = 1.0/delta[0];
        grad_op[(2, 0)] = 1.0/delta[1]; grad_op[(2, 3)] = -1.0/delta[1];
        grad_op[(3, 0)] = -1.0/delta[1]; grad_op[(3, 4)] = 1.0/delta[1];
        grad_op[(4, 0)] = 1.0/delta[2]; grad_op[(4, 5)] = -1.0/delta[2];
        grad_op[(5, 0)] = -1.0/delta[2]; grad_op[(5, 6)] = 1.0/delta[2];
        grad_op
    }
}

