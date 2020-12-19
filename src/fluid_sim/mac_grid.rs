use super::particle::{Particle, Particles};
use super::utils::{sum_from_triplets, Triplets, Grid};
use nalgebra::base::dimension::{U6, U7};
use nalgebra as na;


pub struct MACGrid {
    pub grid: Grid,
    pub density: f32,
    pub q: na::DVector<f32>,
    pub p: na::DVector<f32>,
}

enum Component {
    U,
    V,
    W
}

impl MACGrid {
    pub fn new(grid: &Grid, density: f32) -> Self {
        let Grid {dims, grid_shape} = grid;
        // In order to compute the divergence of the velocity we need information
        // from the current and the three adjacent cells so for cell (i,j,k) we need info 
        // from cell (i+1.j,k), (i.j+1,k) and (i.j,k+1) thus to make the boundary of the grid
        // work out we need.
        let num_velocity_cells = (grid_shape.0 + 1)*(grid_shape.1 + 1)*(grid_shape.2 + 1);
        // On the other hand for pressure we need information for all adjacent cells to compute
        // div grad pressure thus we need a larger grid to handle the boundary nicely.
        let num_pressure_cells = (grid_shape.0 + 2)*(grid_shape.1 + 2)*(grid_shape.2 + 2);
        MACGrid {
            grid: *grid,
            q: na::DVector::<f32>::zeros(num_velocity_cells*3),
            p: na::DVector::<f32>::zeros(num_pressure_cells),
            density
        }
    }

    // Public interface
    pub fn get_velocity(&self, pos: na::Vector3<f32>) -> na::Vector3<f32> {
        unimplemented!();
    }

    pub fn set_from_particles<T>(&mut self, particles: &Particles){
        unimplemented!();
    }

    pub fn pressure_project(&mut self) {
        unimplemented!();
    }

    //Private methods
    // Useful latter for implementing constant flow boundary condition.
    fn boundary_velocity(&self, c:Component, i: usize, j:usize, k: usize) -> f32 {
        0.0
    }

    /// Returns true if the c component of velocity at position i, j, k. 
    /// is outside or on the boundary
    fn is_boundary_velocity(&self,  c:Component, i: usize, j:usize, k: usize) -> bool{
        unimplemented!();
    }

    /// Returns true pressure at point i, j, k is outside or on the boundary.
    fn is_boundary_pressure(&self, i: usize, j:usize, k: usize) -> bool {
        unimplemented!();
    }

    /// Get the the index in q for component_ijk.
    fn get_q_ind(&self, c:Component, i:usize, j:usize, k:usize) -> usize {
        unimplemented!();
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

