pub type Dimensions = [u64; 3]; // TODO: Make into usize? Rather awkward in places right now.
pub type Position = Dimensions;

// TODO: Bitmap optimization.
#[derive(Clone)]
pub struct Mask {
    dimensions: Dimensions,
    cells: Box<[bool]>,
}

impl Mask {
    pub fn new(dimensions: Dimensions) -> Self {
        let [w, h, d] = dimensions.map(|v| v as usize);
        let cells = vec![false; w * h * d].into_boxed_slice();
        Self { dimensions, cells }
    }

    pub fn from_cells(dimensions: [u64; 3], cells: Box<[bool]>) -> Mask {
        assert_eq!(dimensions.iter().product::<u64>(), cells.len() as u64);
        Self { dimensions, cells }
    }

    pub const fn dimensions(&self) -> Dimensions {
        self.dimensions
    }

    pub const fn dimensions_usize(&self) -> [usize; 3] {
        let [x, y, z] = self.dimensions;
        [x as usize, y as usize, z as usize]
    }

    pub fn count<const VALUE: bool>(&self) -> usize {
        self.cells.iter().map(|&c| c as u8 as usize).sum()
    }

    // TODO: Impl get and get_mut as Deref<Slice>? Would that be nicer?
    #[inline(always)]
    pub fn get(&self, idx: Position) -> Option<&bool> {
        let [w, h, _d] = self.dimensions_usize();
        let [x, y, z] = idx.map(|v| v as usize);
        self.cells.get(x + y * w + z * w * h)
    }

    #[inline(always)]
    pub fn get_mut(&mut self, idx: Position) -> Option<&mut bool> {
        let [w, h, _d] = self.dimensions_usize();
        let [x, y, z] = idx.map(|v| v as usize);
        self.cells.get_mut(x + y * w + z * w * h)
    }

    /// Return a three-dimensional index into this [`Mask`] from a linear index.
    ///
    /// This function makes no guarantees about whether the returned index is actually within the
    /// `Mask`.
    #[inline(always)]
    fn idx(&self, mut i: usize) -> Position {
        let [w, h, _d] = self.dimensions_usize();

        let x = i % w;
        i /= w;
        let y = i % h;
        let z = i / h;

        [x as u64, y as u64, z as u64]
    }

    /// Return an [`Iterator`] over all indices where the the cell is equal to `VALUE`.
    pub fn indices_where<const VALUE: bool>(&self) -> impl Iterator<Item = Position> + '_ {
        self.linear_indices_where::<VALUE>().map(|i| self.idx(i))
    }

    /// Return an [`Iterator`] over all linear indices into the internal `cells` where the the cell
    /// is equal to `VALUE`.
    pub fn linear_indices_where<const VALUE: bool>(&self) -> impl Iterator<Item = usize> + '_ {
        self.cells
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v == VALUE { Some(i) } else { None })
    }

    /// Apply some `mask` onto this [`Mask`].
    ///
    /// This is essentially an `or-assign` operation between the `self` and the provided `mask`.
    pub fn apply_mask(&mut self, mask: &Mask) {
        assert_eq!(self.dimensions, mask.dimensions);
        assert_eq!(self.cells.len(), mask.cells.len()); // For good measure, so the compiler gets it.
        self.cells
            .iter_mut()
            .zip(mask.cells.iter())
            .for_each(|(s, &m)| *s |= m);
    }

    pub fn indices_where_into_vec<const VALUE: bool>(&self, locations: &mut Vec<Position>) {
        let [w, h, d] = self.dimensions();
        assert_eq!((w * h * d) as usize, self.cells.len());
        for z in 0..d {
            for y in 0..h {
                for x in 0..w {
                    let lin_idx = x + y * w + z * w * h;
                    // Safety: We know that our access is in bounds since `w*h*d` is equal to
                    // `self.cells.len()`.
                    if unsafe { *self.cells.get_unchecked(lin_idx as usize) } == VALUE {
                        locations.push([x, y, z]);
                    }
                }
            }
        }
    }
}
