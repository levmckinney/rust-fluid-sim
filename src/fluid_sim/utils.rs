use nalgebra as na;

#[derive(Copy, Clone)]
pub struct Grid {
    pub offset: na::Vector3<f32>,
    pub dims: na::Vector3<f32>,
    pub grid_shape: (usize, usize, usize), // Shape of the grid
}

impl Grid {
    pub fn new(offset: na::Vector3<f32>, cell_size: na::Vector3<f32>, grid_shape: (usize, usize, usize)) -> Self {
        Grid {
            offset,
            grid_shape,
            dims: na::Vector3::new(
                cell_size[0]*(grid_shape.0 as f32),
                cell_size[1]*(grid_shape.1 as f32),
                cell_size[2]*(grid_shape.2 as f32))
        }
    }

    pub fn cells(&self) -> itertools::ConsTuples<itertools::Product<itertools::Product<std::ops::Range<usize>, std::ops::Range<usize>>, std::ops::Range<usize>>, ((usize, usize), usize)> {
        let grid_shape = self.grid_shape;
        iproduct!(0..grid_shape.0, 0..grid_shape.1, 0..grid_shape.2)
    }

    pub fn cell_size(&self) -> na::Vector3::<f32> {
        na::Vector3::<f32>::new(
            self.dims[0]/(self.grid_shape.0 as f32), 
            self.dims[1]/(self.grid_shape.1 as f32), 
            self.dims[2]/(self.grid_shape.2 as f32))
    }

    pub fn clamp_inds(&self, i: usize, j:usize, k: usize) -> (usize, usize, usize) {
        return (na::min(i, self.grid_shape.0 - 1),
                na::min(j, self.grid_shape.1 - 1),
                na::min(k, self.grid_shape.2 - 1))
    }

    /// Clamp a position to be with the grid.
    pub fn clamp(&self, position: &na::Vector3<f32>) -> na::Vector3<f32> {
        let epsilon = self.cell_size()*0.01;
        na::Vector3::new(
            na::clamp(position[0], self.offset[0], self.offset[0] + self.dims[0] - epsilon[0]),
            na::clamp(position[1], self.offset[1], self.offset[1] + self.dims[1] - epsilon[1]),
            na::clamp(position[2], self.offset[2], self.offset[2] + self.dims[2] - epsilon[2]))
    }
    
    /// Get the cell at a position within the grid.
    /// Positions outside the grid will be clamped to be within the grid first.
    pub fn get_containing_cell_ind(&self, position: &na::Vector3<f32>) -> (usize, usize, usize) {
        let cell_size = self.cell_size();
        let clamped_centered = self.clamp(position) - self.offset;
        let (i,  j,  k) = (
            (clamped_centered[0] / cell_size[0]).floor() as usize,
            (clamped_centered[1] / cell_size[1]).floor() as usize,
            (clamped_centered[2] / cell_size[2]).floor() as usize);
        assert!(
            i < self.grid_shape.0 && j < self.grid_shape.1 && k < self.grid_shape.2,
            "grid_shape = {:?}, (i, j, k) = {:?}, clamped_centered = {}", 
            self.grid_shape, (i, j, k), clamped_centered
        ); // I messing with a tricky type conversion here.
        (i, j, k)
    }

    /// Get position
    pub fn get_position(&self, i: usize, j: usize, k: usize) -> na::Vector3<f32> {
        let cell_size = self.cell_size();
        let centered = na::Vector3::<f32>::new(
            (i as f32)*cell_size[0], 
            (j as f32)*cell_size[1], 
            (k as f32)*cell_size[2]);
        self.offset + centered
    }

    /// Return wether (i, j, k) is an index of a cell within the grid.
    pub fn index_within(&self, i: usize, j: usize, k: usize) -> bool {
        i < self.grid_shape.0 && j < self.grid_shape.1 && k < self.grid_shape.2
    }

    /// Get number of cells within the grid
    pub fn num_cells(&self) -> usize {
        self.grid_shape.0*self.grid_shape.1*self.grid_shape.2
    }

    /// Get the index in a 1d array corresponding to the 
    /// i, j, k cell in the grid.
    pub fn flat_ind(&self, i: usize, j: usize, k: usize) -> usize {
        let ind = i + j*self.grid_shape.0 + k*self.grid_shape.0*self.grid_shape.1;
        debug_assert!(ind < self.num_cells(), "index out of range");
        ind
    }

    pub fn bilinear_weight(&self, cell_ind: &(usize, usize, usize), position: &na::Vector3<f32>) -> f32 {
        let cell_size = self.cell_size();
        let pos_ind = self.get_containing_cell_ind(position);
        let diff = position - self.get_position(pos_ind.0, pos_ind.1, pos_ind.2) ;
        let dx = diff[0]/cell_size[0];
        let dy = diff[1]/cell_size[1];
        let dz = diff[2]/cell_size[2];
        debug_assert!(!(cell_ind.0  < pos_ind.0 || cell_ind.1  < pos_ind.1 || cell_ind.2  < pos_ind.2), 
          "position not adjacent to cell_ind! cell_ind = {:?} pos_ind = {:?} grid_shape={:?} position={:?}",
          cell_ind, pos_ind, self.grid_shape, position
        );
        let term1 = match cell_ind.0 - pos_ind.0 {
            0 => 1.0 - dx,
            1 => dx,
            _ => panic!("position not adjacent to cell_ind!")
        };
        let term2 = match cell_ind.1 - pos_ind.1 {
            0 => 1.0 - dy,
            1 => dy,
            _ => panic!("position not adjacent to cell_ind!")
        };
        let term3 = match cell_ind.2 - pos_ind.2 {
            0 => 1.0 - dz,
            1 => dz,
            _ => panic!("position not adjacent to cell_ind!")
        };
        term1*term2*term3
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
    for ((i, j), _) in &deduped_triplets {
        assert!(*i < nrows && *j < ncols, "i: {:?}, j:{:?}, nrows: {:?} ncols: {:?}", i, j, nrows, ncols);
    }
    let (inds, vals): (Vec<(usize, usize)>, Vec<f32>) = deduped_triplets.into_iter().unzip();
    let (irows, icols): (Vec<usize>, Vec<usize>) = inds.into_iter().unzip();
    na::CsMatrix::<f32>::from_triplet(nrows, ncols, &irows[..], &icols[..], &vals[..])
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