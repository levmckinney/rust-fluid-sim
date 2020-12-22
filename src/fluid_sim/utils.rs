use nalgebra as na;

#[derive(Copy, Clone)]
pub struct Grid {
    pub offset: na::Vector3<f64>,
    pub dims: na::Vector3<f64>,
    pub grid_shape: (usize, usize, usize), // Shape of the grid
}

impl Grid {
    pub fn new(offset: na::Vector3<f64>, cell_size: na::Vector3<f64>, grid_shape: (usize, usize, usize)) -> Self {
        Grid {
            offset,
            grid_shape,
            dims: na::Vector3::new(
                cell_size[0]*(grid_shape.0 as f64),
                cell_size[1]*(grid_shape.1 as f64),
                cell_size[2]*(grid_shape.2 as f64))
        }
    }

    pub fn cells(&self) -> itertools::ConsTuples<itertools::Product<itertools::Product<std::ops::Range<usize>, std::ops::Range<usize>>, std::ops::Range<usize>>, ((usize, usize), usize)> {
        let grid_shape = self.grid_shape;
        iproduct!(0..grid_shape.0, 0..grid_shape.1, 0..grid_shape.2)
    }

    pub fn cell_size(&self) -> na::Vector3::<f64> {
        na::Vector3::<f64>::new(
            self.dims[0]/(self.grid_shape.0 as f64), 
            self.dims[1]/(self.grid_shape.1 as f64), 
            self.dims[2]/(self.grid_shape.2 as f64))
    }

    pub fn clamp_inds(&self, i: usize, j:usize, k: usize) -> (usize, usize, usize) {
        return (na::min(i, self.grid_shape.0 - 1),
                na::min(j, self.grid_shape.1 - 1),
                na::min(k, self.grid_shape.2 - 1))
    }

    /// Clamp a position to be with the grid.
    pub fn clamp(&self, position: &na::Vector3<f64>) -> na::Vector3<f64> {
        let epsilon = self.cell_size()*0.01;
        na::Vector3::new(
            na::clamp(position[0], self.offset[0], self.offset[0] + self.dims[0] - epsilon[0]),
            na::clamp(position[1], self.offset[1], self.offset[1] + self.dims[1] - epsilon[1]),
            na::clamp(position[2], self.offset[2], self.offset[2] + self.dims[2] - epsilon[2]))
    }
    
    /// Get the cell at a position within the grid.
    /// Positions outside the grid will be clamped to be within the grid first.
    pub fn get_containing_cell_ind(&self, position: &na::Vector3<f64>) -> (usize, usize, usize) {
        let cell_size = self.cell_size();
        let clamped_centered = self.clamp(position) - self.offset;
        let (i,  j,  k) = (
            (clamped_centered[0] / cell_size[0]).floor() as usize,
            (clamped_centered[1] / cell_size[1]).floor() as usize,
            (clamped_centered[2] / cell_size[2]).floor() as usize);
        assert!(
            i < self.grid_shape.0 && j < self.grid_shape.1 && k < self.grid_shape.2,
            "grid_shape = {:?}, (i, j, k) = {:?}, clamped_centered = {}", 
            self.grid_shape, (i, j, k), clamped_centered); // I messing with a tricky type conversion here.
        (i, j, k)
    }

    /// Get position
    pub fn get_position(&self, i: usize, j: usize, k: usize) -> na::Vector3<f64> {
        let cell_size = self.cell_size();
        let centered = na::Vector3::<f64>::new(
            (i as f64)*cell_size[0], 
            (j as f64)*cell_size[1], 
            (k as f64)*cell_size[2]);
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

    pub fn bilinear_weight(&self, cell_ind: &(usize, usize, usize), position: &na::Vector3<f64>) -> f64 {
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
