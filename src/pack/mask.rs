pub type Dimensions = [u64; 3]; // TODO: Make into usize? Rather awkward in places right now.
pub type Position = Dimensions;

type Backing = u8;
const BACKING_BITS: usize = Backing::BITS as usize;

#[derive(Debug, Clone)]
pub struct Mask {
    dimensions: [usize; 3],
    backings: Box<[Backing]>,
}

impl Mask {
    pub fn new(dimensions: Dimensions) -> Self {
        let dimensions = dimensions.map(|v| v as usize);
        let n: usize = dimensions.iter().product();
        let n_backings = n.div_ceil(BACKING_BITS);
        Self {
            dimensions,
            backings: vec![0; n_backings].into_boxed_slice(),
        }
    }

    /// # Panics
    ///
    /// If the provided `dimensions` are not compatible with the number of `cells`, this function
    /// will panic.
    pub fn from_cells(dimensions: Dimensions, cells: &[bool]) -> Self {
        let mut mask = Self::new(dimensions);
        assert_eq!(
            mask.n_cells(),
            cells.len(),
            "the number of cells must match the dimensions"
        );

        for lin_idx in cells
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v { Some(i) } else { None })
        {
            mask.set_linear_unchecked::<true>(lin_idx)
        }

        mask
    }

    pub const fn dimensions(&self) -> Dimensions {
        let [w, h, d] = self.dimensions;
        [w as u64, h as u64, d as u64]
    }

    /// Return the number of cells that are equal to `VALUE`.
    pub fn count<const VALUE: bool>(&self) -> usize {
        // If we're counting the `false` cells, we may overcount when there are trailing bits in
        // the last backing.
        let overshoot = if VALUE {
            0
        } else {
            self.n_backings() * BACKING_BITS - self.n_cells()
        };
        let uncorrected: usize = self
            .backings
            .iter()
            .map(|backing| {
                if VALUE {
                    backing.count_ones() as usize
                } else {
                    backing.count_zeros() as usize
                }
            })
            .sum();
        uncorrected - overshoot
    }

    const fn n_cells(&self) -> usize {
        let [w, h, d] = self.dimensions;
        w * h * d
    }

    const fn n_backings(&self) -> usize {
        self.backings.len()
    }

    const fn backing_idx(&self, lin_idx: usize) -> (usize, usize) {
        (lin_idx / BACKING_BITS, lin_idx % BACKING_BITS)
    }

    const fn contains(&self, idx: Position) -> bool {
        let [x, y, z] = idx;
        let [w, h, d] = self.dimensions();
        x < w && y < h && z < d
    }

    const fn linear_idx(&self, idx: Position) -> usize {
        let [x, y, z] = idx;
        let [w, h, _] = self.dimensions;
        let lin_idx = x as usize + y as usize * w + z as usize * w * h;
        debug_assert!(lin_idx < self.n_cells());
        lin_idx
    }

    const fn get_linear_unchecked(&self, lin_idx: usize) -> bool {
        debug_assert!(lin_idx < self.n_cells());
        let (backing_idx, bit_idx) = self.backing_idx(lin_idx);
        self.backings[backing_idx] >> bit_idx & 1 != 0
    }

    /// Check whether the cell at the specified `lin_idx` is `VALUE`.
    const fn query_linear_unchecked<const VALUE: bool>(&self, lin_idx: usize) -> bool {
        debug_assert!(lin_idx < self.n_cells());
        let (backing_idx, bit_idx) = self.backing_idx(lin_idx);
        let backing = self.backings[backing_idx];
        // FIXME: This seems convoluted for what the truth table actually looks like.
        if VALUE {
            match backing {
                0 => false,
                u8::MAX => true,
                _ => backing >> bit_idx & 1 != 0,
            }
        } else {
            match backing {
                0 => true,
                u8::MAX => false,
                _ => backing >> bit_idx & 1 == 0,
            }
        }
    }

    fn set_linear_unchecked<const VALUE: bool>(&mut self, lin_idx: usize) {
        debug_assert!(lin_idx < self.n_cells());
        let (backing_idx, bit_idx) = self.backing_idx(lin_idx);
        if VALUE {
            let bit = 1 << bit_idx;
            self.backings[backing_idx] |= bit;
        } else {
            let backing = &mut self.backings[backing_idx];
            *backing &= !(1 << bit_idx);
        }
    }

    pub const fn get(&self, idx: Position) -> Option<bool> {
        if self.contains(idx) {
            let lin_idx = self.linear_idx(idx);
            Some(self.get_linear_unchecked(lin_idx))
        } else {
            None
        }
    }

    /// # Panics
    ///
    /// If `idx` is not within the dimensions of this `Mask`, this function will panic.
    pub fn set(&mut self, idx: Position, value: bool) {
        assert!(
            self.contains(idx),
            "the provided `idx` must be within the dimensions of this mask"
        );
        let lin_idx = self.linear_idx(idx);
        if value {
            self.set_linear_unchecked::<true>(lin_idx)
        } else {
            self.set_linear_unchecked::<false>(lin_idx)
        }
    }

    /// Return an [`Iterator`] over all indices where the cell is equal to `VALUE`.
    pub fn indices_where<const VALUE: bool>(&self) -> impl Iterator<Item = Position> + '_ {
        let [w, h, d] = self.dimensions();
        (0..d).flat_map(move |z| {
            (0..h).flat_map(move |y| {
                (0..w).filter_map(move |x| {
                    let idx = [x, y, z];
                    let lin_idx = self.linear_idx(idx);
                    if self.query_linear_unchecked::<VALUE>(lin_idx) {
                        Some(idx)
                    } else {
                        None
                    }
                })
            })
        })
    }

    /// Apply some `mask` onto this [`Mask`].
    ///
    /// This is essentially an `or-assign` operation between the `self` and the provided `mask`.
    /// # Panics
    ///
    /// If the dimensions of the stamped and stamping mask are different, this function will panic.
    pub fn apply_mask(&mut self, mask: &Mask) {
        assert_eq!(
            self.dimensions(),
            mask.dimensions(),
            "the dimensions of both masks must be identical"
        );
        // For good measure, so the compiler gets it.
        assert_eq!(self.n_backings(), mask.n_backings()); // FIXME: Is this one necessary?

        self.backings
            .iter_mut()
            .zip(mask.backings.iter())
            .for_each(|(s, &m)| *s |= m);
    }
}
