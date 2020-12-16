use nalgebra as na;

/// Create sparse matrix from triplets. 
/// Sum duplicate entries.
fn sum_from_triplets<T: na::Scalar>(triplets: Vec<((usize, usize), T)>) -> na::CsMatrix<T> {
    unimplemented!();
}