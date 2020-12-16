use super::particle::{Particle};
use nalgebra as na;

pub struct MACGrid {
    // Generalized velocity
    dims: na::Vector3<f32>,
    grid_shape: (usize, usize, usize),
    density: f32,
    div_grad_op: na::CsMatrix<f32>,
    div_op: na::CsMatrix<f32>,
    q: na::DVector<f32>,
    p: na::DVector<f32>,
}

enum Component {
    U,
    V,
    W
}

impl MACGrid {
    pub fn new(dims: na::Vector3<f32>, grid_shape: (usize, usize, usize), density: f32) -> Self {
        let num_cells = grid_shape.0*grid_shape.1*grid_shape.2;
        // We need to add one here because we must have 
        let num_velocities = (grid_shape.0 + 1)*(grid_shape.1 + 1)*(grid_shape.2 + 1);
        let deltas = na::Vector3::<f32>::new(dims[0]/(grid_shape.0 as f32), 
                                           dims[1]/(grid_shape.1 as f32), 
                                           dims[2]/(grid_shape.2 as f32))
        MACGrid{ 
            dims,
            grid_shape,
            div_grad_op: div_of_grad_operator_global(grid_shape, deltas),
            div_op: div_operator_global(grid_shape, deltas),
            density,
            q: na::DVector::<f32>::zeros(num_velocities*3),
            p: na::DVector::<f32>::zeros(num_cells)
        }
    }

    // Public interface
    pub fn get_velocity(&self, pos: na::Vector3<f32>) -> na::Vector3<f32> {
        unimplemented!();
    }

    pub fn set_from_particles<T>(&mut self, particles: &[Particle<T>]){
        unimplemented!();
    }

    pub fn pressure_project(&mut self) {
        unimplemented!();
    }

    //Private methods
    /// Returns the pressure at this point on the boundary get pressure at boundary
    fn boundary_condition_pressure(&self, i: usize, j:usize, k: usize) -> f32{
        0.0
    }

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
    fn get_boundary_filter(&self, i:usize, j:usize, k:usize) -> na::MatrixN<f32,6> {
        unimplemented!();
    }

    /// Returns local div operator
    fn div_operator_local(&self, deltas: na::Vector3<f32>) -> na::Matrix3x1<f32> {
        unimplemented!();
    }

    // what is the shape of the input operator vector?
    fn div_of_grad_operator_local(&self, deltas: na::Vector3<f32>) -> na::MatrixN<f32,6> {
        unimplemented!();
    }

    /// Returns global divergence operator matrix.
    fn div_operator_global(&self, deltas: na::Vector3<f32>) -> na::CsMatrix<f32>{
        unimplemented!();
    } 

    /// Returns global divergence grad operator.
    fn div_of_grad_operator_global(&self, deltas: na::Vector3<f32>) -> na::CsMatrix<f32>{
        unimplemented!();
    }
}

