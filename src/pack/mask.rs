pub type Dimensions = [u64; 3]; // TODO: Make into usize? Rather awkward in places right now.
pub type Position = Dimensions;
pub type SignedPosition = [i64; 3];

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
        assert!(
            dimensions.iter().all(|&v| v > 0),
            "a mask must have non-zero dimensions"
        );
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
        assert_eq!(
            dimensions.iter().product::<u64>() as usize,
            cells.len(),
            "the number of cells must match the dimensions"
        );

        Self::from_cells_iter(dimensions, cells.into_iter().copied())
    }

    /// An internal function for setting up a [`Mask`] from an iterator of `bool`s.
    ///
    /// Providing an iterator that is not exactly the correct length is not acceptable behavior but
    /// is not unsafe.
    fn from_cells_iter(dimensions: Dimensions, cells: impl Iterator<Item = bool>) -> Self {
        let mut mask = Self::new(dimensions);

        for lin_idx in cells
            .enumerate()
            .filter_map(|(i, v)| if v { Some(i) } else { None })
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

    pub const fn spatial_idx(&self, mut lin_idx: usize) -> Option<Position> {
        if lin_idx >= self.n_cells() {
            return None;
        }

        let [w, h, _] = self.dimensions;
        let x = lin_idx % w;
        lin_idx /= w;
        let y = lin_idx % h;
        lin_idx /= h;
        let z = lin_idx;

        Some([x as u64, y as u64, z as u64])
    }

    const fn linear_idx(&self, idx: Position) -> usize {
        let [x, y, z] = idx;
        let [w, h, _] = self.dimensions;
        let lin_idx = x as usize + y as usize * w + z as usize * w * h;
        debug_assert!(lin_idx < self.n_cells());
        lin_idx
    }

    const fn linear_idx_periodic(&self, idx: Position) -> usize {
        self.linear_idx(normalize_periodic(idx, self.dimensions()))
    }

    const fn linear_idx_periodic_signed(&self, idx: SignedPosition) -> usize {
        self.linear_idx(normalize_periodic_signed(idx, self.dimensions()))
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
                Backing::MAX => true,
                _ => backing >> bit_idx & 1 != 0,
            }
        } else {
            match backing {
                0 => true,
                Backing::MAX => false,
                _ => backing >> bit_idx & 1 == 0,
            }
        }
    }

    fn set_linear_unchecked<const VALUE: bool>(&mut self, lin_idx: usize) {
        debug_assert!(lin_idx < self.n_cells());
        let (backing_idx, bit_idx) = self.backing_idx(lin_idx);
        if VALUE {
            self.backings[backing_idx] |= 1 << bit_idx;
        } else {
            self.backings[backing_idx] &= !(1 << bit_idx);
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

    pub const fn get_periodic(&self, idx: Position) -> bool {
        let lin_idx = self.linear_idx(normalize_periodic(idx, self.dimensions()));
        self.get_linear_unchecked(lin_idx)
    }

    pub const fn get_periodic_signed(&self, idx: SignedPosition) -> bool {
        let lin_idx = self.linear_idx(normalize_periodic_signed(idx, self.dimensions()));
        self.get_linear_unchecked(lin_idx)
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

    pub fn set_periodic(&mut self, idx: Position, value: bool) {
        let lin_idx = self.linear_idx(normalize_periodic(idx, self.dimensions()));
        if value {
            self.set_linear_unchecked::<true>(lin_idx)
        } else {
            self.set_linear_unchecked::<false>(lin_idx)
        }
    }

    pub fn set_periodic_signed(&mut self, idx: SignedPosition, value: bool) {
        let lin_idx = self.linear_idx(normalize_periodic_signed(idx, self.dimensions()));
        if value {
            self.set_linear_unchecked::<true>(lin_idx)
        } else {
            self.set_linear_unchecked::<false>(lin_idx)
        }
    }

    /// Return an [`Iterator`] over all linear indices where the cell is equal to `VALUE`.
    pub fn linear_indices_where<const VALUE: bool>(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.n_cells()).filter(|&lin_idx| self.query_linear_unchecked::<VALUE>(lin_idx))
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

/// Take an `idx` that may fall outside of the `dimensions` and return a normalized [`Position`].
const fn normalize_periodic(idx: Position, dimensions: Dimensions) -> Position {
    // FIXME: This const unrolling sucks. See if I can do something about that.
    [
        idx[0] % dimensions[0],
        idx[1] % dimensions[1],
        idx[2] % dimensions[2],
    ]
}

/// Take a signed `idx` that may fall outside of the `dimensions` and return a normalized [`Position`].
const fn normalize_periodic_signed(idx: SignedPosition, dimensions: Dimensions) -> Position {
    // FIXME: This const unrolling sucks. See if I can do something about that.
    [
        idx[0].rem_euclid(dimensions[0] as i64) as u64,
        idx[1].rem_euclid(dimensions[1] as i64) as u64,
        idx[2].rem_euclid(dimensions[2] as i64) as u64,
    ]
}

pub struct MaskSlice<'m> {
    start: Position,
    dimensions: Dimensions,
    mask: &'m Mask,
}

impl<'m> Mask {
    pub fn slice_from_start_dimensions(
        &'m self,
        start: Position,
        dimensions: Dimensions,
    ) -> MaskSlice<'m> {
        // TODO: A bit sloppy all of this math. Maybe we can just calculate these things with a U64Vec3?
        let [w, h, d] = dimensions;
        let [sx, sy, sz] = start;
        let end = [sx + w, sy + h, sz + d];
        {
            let [w, h, d] = self.dimensions(); // The dimensions of the underlying Mask.
            let [sx, sy, sz] = start;
            assert!(sx < w && sy < h && sz < d);
            let [ex, ey, ez] = end;
            assert!(ex < w && ey < h && ez < d);
        }
        assert!(dimensions.iter().all(|&v| v > 0));

        MaskSlice {
            start,
            dimensions,
            mask: &self,
        }
    }

    pub fn slice(&'m self, ranges: [std::ops::Range<u64>; 3]) -> MaskSlice<'m> {
        let [w, h, d] = self.dimensions(); // The dimensions of the underlying Mask.

        let start = ranges.clone().map(|r| r.start);
        let end = ranges.clone().map(|r| r.end);
        {
            let [sx, sy, sz] = start;
            assert!(sx < w && sy < h && sz < d);
            let [ex, ey, ez] = end;
            assert!(ex < w && ey < h && ez < d);
        }
        let dimensions = ranges.map(|r| r.end - r.start);
        assert!(dimensions.iter().all(|&v| v > 0));

        MaskSlice {
            start,
            dimensions,
            mask: &self,
        }
    }
}

impl MaskSlice<'_> {
    pub fn iter_linear(&self) -> impl Iterator<Item = bool> + '_ {
        let [start_x, start_y, start_z] = self.start;
        let [w, h, d] = self.dimensions;
        let [end_x, end_y, end_z] = [start_x + w, start_y + h, start_z + d];

        (start_z..end_z).flat_map(move |z| {
            (start_y..end_y).flat_map(move |y| {
                (start_x..end_x).map(move |x| {
                    // FIXME: We know that this is valid. Time for a get_unchecked?
                    self.mask.get([x, y, z]).unwrap()
                })
            })
        })
    }

    /// Iterator over the [`Position`]s of cells equal to `VALUE` in this [`MaskSlice`].
    ///
    /// The `Positions` are according to the internal `Mask`'s absolute coordinates.
    pub fn iter_where<const VALUE: bool>(&self) -> impl Iterator<Item = Position> + '_ {
        let [start_x, start_y, start_z] = self.start;
        let [w, h, d] = self.dimensions;
        let [end_x, end_y, end_z] = [start_x + w, start_y + h, start_z + d];

        (start_z..end_z).flat_map(move |z| {
            (start_y..end_y).flat_map(move |y| {
                (start_x..end_x).flat_map(move |x| {
                    let idx = [x, y, z];
                    let lin_idx = self.mask.linear_idx(idx);
                    if self.mask.query_linear_unchecked::<VALUE>(lin_idx) {
                        Some(idx)
                    } else {
                        None
                    }
                })
            })
        })
    }

    pub fn to_mask(&self) -> Mask {
        let cells = self.iter_linear();
        Mask::from_cells_iter(self.dimensions, cells)
    }
}

impl MaskSlice<'_> {
    // FIXME: Probably not needed.
    fn all<const VALUE: bool>(&self) -> bool {
        self.iter_linear().all(|cell| cell == VALUE)
    }

    pub fn any<const VALUE: bool>(&self) -> bool {
        self.iter_linear().any(|cell| cell == VALUE)
    }

    // FIXME: Probably not needed.
    /// Return the [`Position`] of the first cell equal to `VALUE` in this [`MaskSlice`] if it exists.
    ///
    /// The `Positions` are relative to the `MaskSlice` `start`.
    fn first_where<const VALUE: bool>(&self) -> Option<Position> {
        self.iter_where::<VALUE>().next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create() {
        let _ = Mask::new([10, 20, 30]);
    }

    #[test]
    #[should_panic]
    fn create_zero() {
        let _ = Mask::new([0, 0, 0]);
    }

    #[test]
    fn from_cells() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let _ = Mask::from_cells([3, 3, 3], &cells);
    }

    #[test]
    #[should_panic]
    fn from_cells_under() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true,
        ];
        let _ = Mask::from_cells([3, 3, 3], &cells);
    }

    #[test]
    #[should_panic]
    fn from_cells_over() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true, true, true, false, true,
        ];
        let _ = Mask::from_cells([3, 3, 3], &cells);
    }

    #[test]
    fn dimensions() {
        let dimensions = [5172, 1312, 161];
        let mask = Mask::new(dimensions);
        assert_eq!(mask.dimensions(), dimensions);
    }

    #[test]
    fn count() {
        let mut mask = Mask::new([172, 1312, 161]);

        mask.set([123, 456, 78], true);
        mask.set([123, 784, 56], true);
        mask.set([78, 1236, 45], true);

        assert_eq!(mask.count::<true>(), 3);
        assert_eq!(mask.count::<false>(), mask.n_cells() - mask.count::<true>());

        mask.set([78, 1236, 45], false);
        assert_eq!(mask.count::<true>(), 2);
        assert_eq!(mask.count::<false>(), mask.n_cells() - mask.count::<true>());
    }

    #[test]
    fn n_cells() {
        let mask = Mask::new([123, 456, 789]);
        assert_eq!(mask.n_cells(), 123 * 456 * 789);
    }

    #[test]
    fn n_backings() {
        let mask = Mask::new([123, 456, 789]);
        assert_eq!(mask.n_backings(), mask.n_cells().div_ceil(BACKING_BITS));
    }

    #[test]
    fn backing_idx() {
        let mask = Mask::new([123, 456, 789]);
        assert_eq!(mask.backing_idx(0), (0, 0));
        assert_eq!(mask.backing_idx(1), (0, 1));
        assert_eq!(mask.backing_idx(2), (0, 2));
        assert_eq!(mask.backing_idx(8), (8 / BACKING_BITS, 8 % BACKING_BITS));
    }

    #[test]
    fn contains() {
        let mask = Mask::new([123, 456, 789]);
        assert!(!mask.contains([123, 456, 789]));
        assert!(!mask.contains([122, 456, 789]));
        assert!(!mask.contains([122, 455, 789]));
        assert!(mask.contains([122, 455, 788]));
        assert!(mask.contains([0, 0, 0]));
        assert!(!mask.contains([12345, 0, 0]));
    }

    #[test]
    fn spatial_idx() {
        let mask = Mask::new([100, 100, 100]);
        assert_eq!(mask.spatial_idx(0), Some([0, 0, 0]));
        assert_eq!(mask.spatial_idx(99), Some([99, 0, 0]));
        assert_eq!(mask.spatial_idx(100), Some([0, 1, 0]));
        assert_eq!(mask.spatial_idx(101), Some([1, 1, 0]));
        assert_eq!(mask.spatial_idx(1000000 - 1), Some([99, 99, 99]));
        assert!(mask.spatial_idx(1000000).is_none());
        assert!(mask.spatial_idx(12345678).is_none());
    }

    #[test]
    fn linear_idx() {
        let mask = Mask::new([100, 100, 100]);
        assert_eq!(mask.linear_idx([0, 0, 0]), 0);
        assert_eq!(mask.linear_idx([99, 0, 0]), 99);
        assert_eq!(mask.linear_idx([0, 1, 0]), 100);
        assert_eq!(mask.linear_idx([1, 1, 0]), 101);
        assert_eq!(mask.linear_idx([99, 99, 99]), 1000000 - 1);
    }

    #[test]
    #[should_panic]
    fn linear_idx_panic() {
        let mask = Mask::new([100, 100, 100]);
        assert_eq!(mask.linear_idx([100, 100, 100]), 1000000); // FIXME: Reconsider this behavior.
    }

    #[test]
    fn query_linear_unchecked() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mask = Mask::from_cells([3, 3, 3], &cells);

        assert!(mask.query_linear_unchecked::<false>(0));
        assert!(mask.query_linear_unchecked::<true>(13));
        assert!(mask.query_linear_unchecked::<false>(7));
        assert!(mask.query_linear_unchecked::<true>(12));
    }

    #[test]
    fn set_linear_unchecked() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mut mask = Mask::from_cells([3, 3, 3], &cells);

        // Was false, set false.
        assert!(mask.query_linear_unchecked::<false>(0));
        mask.set_linear_unchecked::<false>(0);
        assert!(mask.query_linear_unchecked::<false>(0));

        // Was true, set false.
        assert!(mask.query_linear_unchecked::<true>(13));
        mask.set_linear_unchecked::<false>(13);
        assert!(mask.query_linear_unchecked::<false>(13));

        // Was false, set true.
        assert!(mask.query_linear_unchecked::<false>(8));
        mask.set_linear_unchecked::<true>(8);
        assert!(mask.query_linear_unchecked::<true>(8));

        // Was true, set true.
        assert!(mask.query_linear_unchecked::<true>(12));
        mask.set_linear_unchecked::<true>(12);
        assert!(mask.query_linear_unchecked::<true>(12));

        // Was true, set false.
        assert!(mask.query_linear_unchecked::<true>(17));
        mask.set_linear_unchecked::<false>(17);
        assert!(mask.query_linear_unchecked::<false>(17));
    }

    #[test]
    fn get() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mut mask = Mask::from_cells([3, 3, 3], &cells);

        assert_eq!(mask.get([0, 0, 0]), Some(false));
        assert_eq!(mask.get([1, 1, 1]), Some(true));
        assert_eq!(mask.get([2, 2, 2]), Some(true));
        assert_eq!(mask.get([2, 2, 3]), None);
        assert_eq!(mask.get([3, 3, 3]), None);
        assert_eq!(mask.get([2, 4, 2]), None);

        mask.set([1, 1, 1], true);
        assert_eq!(mask.get([1, 1, 1]), Some(true));
    }

    #[test]
    fn set() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mut mask = Mask::from_cells([3, 3, 3], &cells);

        // Was false, set false.
        assert_eq!(mask.get([0, 0, 0]), Some(false));
        mask.set([0, 0, 0], false);
        assert_eq!(mask.get([0, 0, 0]), Some(false));

        // Was true, set true.
        assert_eq!(mask.get([0, 1, 0]), Some(true));
        mask.set([0, 1, 0], true);
        assert_eq!(mask.get([0, 1, 0]), Some(true));

        // Was false, set true.
        assert_eq!(mask.get([1, 2, 0]), Some(false));
        mask.set([1, 2, 0], true);
        assert_eq!(mask.get([1, 2, 0]), Some(true));

        // Was true, set false.
        assert_eq!(mask.get([2, 2, 2]), Some(true));
        mask.set([2, 2, 2], false);
        assert_eq!(mask.get([2, 2, 2]), Some(false));
    }

    #[test]
    fn indices_where() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mask = Mask::from_cells([3, 3, 3], &cells);

        let where_false = [
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [0, 2, 0],
            [1, 2, 0],
            [2, 2, 0],
            [0, 0, 1],
            [2, 0, 1],
            [0, 2, 1],
            [1, 2, 1],
            [0, 0, 2],
            [1, 0, 2],
            [2, 0, 2],
        ];

        assert_eq!(
            mask.indices_where::<false>().collect::<Vec<_>>(),
            where_false
        );

        let where_true = [
            [0, 1, 0],
            [1, 1, 0],
            [2, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
            [2, 1, 1],
            [2, 2, 1],
            [0, 1, 2],
            [1, 1, 2],
            [2, 1, 2],
            [0, 2, 2],
            [1, 2, 2],
            [2, 2, 2],
        ];

        assert_eq!(mask.indices_where::<true>().collect::<Vec<_>>(), where_true);
    }

    #[test]
    fn linear_indices_where() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mask = Mask::from_cells([3, 3, 3], &cells);

        let where_false = [0, 1, 2, 6, 7, 8, 9, 11, 15, 16, 18, 19, 20];

        assert_eq!(
            mask.linear_indices_where::<false>().collect::<Vec<_>>(),
            where_false
        );

        let where_true = [3, 4, 5, 10, 12, 13, 14, 17, 21, 22, 23, 24, 25, 26];

        assert_eq!(
            mask.linear_indices_where::<true>().collect::<Vec<_>>(),
            where_true
        );
    }

    #[test]
    fn apply_mask() {
        let cells = [
            false, false, false, true, true, true, false, false, false, false, true, false, true,
            true, true, false, false, true, false, false, false, true, true, true, true, true,
            true,
        ];
        let mut base = Mask::from_cells([3, 3, 3], &cells);
        assert!(base.query_linear_unchecked::<false>(0));
        assert!(base.query_linear_unchecked::<false>(1));

        let mut mask = base.clone();
        mask.set([1, 0, 0], true);
        assert!(mask.query_linear_unchecked::<false>(0));
        assert!(mask.query_linear_unchecked::<true>(1));

        // Mask the base.
        base.apply_mask(&mask);
        // And now base should have that same bit set and be identical aside from that.
        assert!(base.query_linear_unchecked::<false>(0));
        assert!(base.query_linear_unchecked::<true>(1));
    }
}