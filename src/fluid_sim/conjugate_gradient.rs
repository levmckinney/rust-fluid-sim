use super::mac_grid::MACGrid;
use nalgebra as na;

pub fn conjugate_gradient(
    mac_grid: &MACGrid,
    init_guess: &na::DVector<f64>,
    b: &na::DVector<f64>,
    tol: f64,
    max_iterations: usize) -> na::DVector<f64> {
    let mut x = init_guess.clone();
    let mut r = b - mac_grid.div_grad_pressure(init_guess); // TODO put this behind a trait
    if r.abs().max() < tol { // Init guess good enough
        return x;
    }
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
    return x;
}