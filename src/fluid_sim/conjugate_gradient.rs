use super::mac_grid::MACGrid;
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
    mac_grid: &MACGrid,
    init_guess: &na::DVector<f32>,
    b: &na::DVector<f32>,
    tol: f32,
    max_iterations: usize) -> na::DVector<f32> {
    let mut x = init_guess.clone();
    let mut r = b - mac_grid.div_grad_pressure(init_guess); // TODO put this behind a trait
    if r.abs().max() < tol { // Init guess good enough
        return x;
    }
    println!("before max div {}", r.abs().max());
    let mut p = r.clone();
    let mut rsold = r.dot(&r);

    for _ in 0..max_iterations  {
        let sprs_mat_x_p = mac_grid.div_grad_pressure(&p);
        let alpha = rsold/p.dot(&sprs_mat_x_p);
        x += alpha*&p;
        r -= alpha*sprs_mat_x_p;
        if r.abs().max() < tol { // Converged!
            break; 
        }
        let rsnew = r.dot(&r);
        let beta = rsnew/rsold;
        p = &r + beta*p;
        rsold = rsnew;
    }
    println!("max div {}", (b - mac_grid.div_grad_pressure(&x)).abs().max());
    return x;
}