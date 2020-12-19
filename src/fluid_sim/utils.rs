use nalgebra as na;

#[derive(Copy, Clone)]
pub struct Grid {
    pub dims: na::Vector3<f32>,
    pub grid_shape: (usize, usize, usize),
}

impl Grid {
    pub fn cells(&self) -> itertools::ConsTuples<itertools::Product<itertools::Product<std::ops::Range<usize>, std::ops::Range<usize>>, std::ops::Range<usize>>, ((usize, usize), usize)> {
        let grid_shape = self.grid_shape;
        iproduct!(0..grid_shape.0, 0..grid_shape.1, 0..grid_shape.2)
    }

    pub fn cell_size(&self) -> na::Vector3::<f32> {
        let Grid {dims, grid_shape} = self;
        return na::Vector3::<f32>::new(
            dims[0]/(grid_shape.0 as f32), 
            dims[1]/(grid_shape.1 as f32), 
            dims[2]/(grid_shape.2 as f32));
    }

    pub fn within(&self, i: usize, j: usize, k: usize) -> bool {
        return i < self.grid_shape.0 && j < self.grid_shape.1 && k < self.grid_shape.2
    }
}

pub type Triplet = ((usize, usize), f32);
pub type Triplets = Vec<Triplet>;

fn sum_duplicate_triplets(triplets: Triplets) -> Triplets {
    let mut triplets = triplets.clone();
    triplets.sort_by(|item1, item2|
        if item1.0.0 == item2.0.0 {
            item1.0.1.cmp(&item2.0.1)
        } else {
            item1.0.0.cmp(&item2.0.0)
        }
    );

    let mut deduped_triplets = Vec::new();
    let mut last = None;
    let mut sum = 0.0;
    for (indices, value) in triplets {
        if last.is_none() || indices == last.unwrap() {
            sum += value;
        } else if last.is_some(){
            deduped_triplets.push((last.unwrap(), sum));
            sum = value;
        }
        last = Some(indices);
    }
    deduped_triplets.push((last.unwrap(), sum));

    deduped_triplets
}

/// Create sparse matrix from triplets. 
/// Sum duplicate entries.
pub fn sum_from_triplets(nrows: usize, ncols: usize, triplets: Triplets) -> na::CsMatrix<f32> {
    let deduped_triplets = sum_duplicate_triplets(triplets);
    let (inds, vals): (Vec<(usize, usize)>, Vec<f32>) = deduped_triplets.into_iter().unzip();
    let (irows, icols): (Vec<usize>, Vec<usize>) = inds.into_iter().unzip();
    na::CsMatrix::<f32>::from_triplet(nrows, ncols, &irows[..], &icols[..], &vals[..])
}


/// Axis aligned bounding box.
pub struct AABB {
    dims: na::Vector3<f32>,
    corner: na::Vector3<f32>
}

impl AABB {
    /// Compute the trilinear weight within the box.
    pub fn trilinear_weight(position: na::Vector3<f32>) -> f32 {
        unimplemented!();
    }
}

pub fn clamp(value: na::Vector3<f32>, lower: na::Vector3<f32>, upper: na::Vector3<f32>) -> na::Vector3<f32> {
    na::Vector3::new(
        na::clamp(value[0], lower[0], upper[0]),
        na::clamp(value[1], lower[1], upper[1]),
        na::clamp(value[2], lower[2], upper[2]))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sum_duplicate_triplets() {
        assert_eq!(sum_duplicate_triplets(
            vec![((1, 2), 2.0),
            ((1, 1), 0.6),
            ((1, 1), 0.4)]),
            vec![((1, 1), 1.0), ((1, 2), 2.0)]           
        )
    }
}