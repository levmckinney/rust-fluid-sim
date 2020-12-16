use nalgebra as na;

/// Create the precondioner matrix for a fluid simulation we use the same one as in
/// the fluid notes. #TODO add name of it here
pub fn fluid_sim_preconditioner() -> na::CsMatrix{
    unimplemented!();
}

/// Run preconditioned conjugate gradient algorithm decomposition
pub fn pre_conditioned_conjugate_gradient(
    pre_cond: na::CsMatrix, 
    A: na::CsMatrix,
    y: na::VectorD,
    iterations: usize,
) -> na::VectorD<f32> {
    unimplemented!();
}