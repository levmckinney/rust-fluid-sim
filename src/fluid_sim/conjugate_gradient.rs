use nalgebra as na;

/// Create the precondioner matrix for a fluid simulation we use the same one as in
/// the fluid notes. #TODO add name of it here
pub fn fluid_sim_preconditioner() -> na::CsMatrix<f32> {
    unimplemented!();
}

/// Run preconditioned conjugate gradient algorithm decomposition
pub fn pre_conditioned_conjugate_gradient(
    pre_cond: na::CsMatrix<f32>, 
    A: na::CsMatrix<f32>,
    y: na::DVector<f32>,
    iterations: usize,
) -> na::DVector<f32> {
    unimplemented!();
}

pub fn conjugate_gradient(
    sprs_mat: &na::CsMatrix<f32>,
    init_guess: &na::DVector<f32>,
    b: &na::DVector<f32>,
    tol: f32,
    max_iterations: usize) -> na::DVector<f32> {
    let mut x = init_guess;
    let mut r = b - sprs_mat*x.into();
    let mut p = r;
    let mut i = 0.0;

    while i < max_iterations && r.abs().max() > tol {
        let sprs_mat_x_p = sprs_mat*p;
        let alpha = r.dot(r)/(p.transpose()*sprs_mat_x_p);
        x += alpha*p;
        r -= alpha*sprs_mat_x_p;
    }
}