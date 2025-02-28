use eightyseven::structure::Atom;
use eightyseven::writer::WriteGro;
use glam::{UVec3, Vec3};

use crate::{
    convert::Convert,
    structure::{BoxVecsExtension, Structure},
};

/// Return an index into a linear array that represents items on a 3-dimensional grid in a z-major
/// ordering.
pub fn index_3d(pos: UVec3, dimensions: UVec3) -> usize {
    // Casts to usize here serve to avoid wrapping when indexing _very_ large sizes. You never know!
    pos.x as usize
        + pos.y as usize * dimensions.x as usize
        + pos.z as usize * dimensions.x as usize * dimensions.y as usize
}

/// Iterate over all positions between `(0, 0, 0)` and `dimensions`.
pub fn iter_3d(dimensions: UVec3) -> impl Iterator<Item = UVec3> {
    let [dx, dy, dz] = dimensions.to_array();
    (0..dz).flat_map(move |z| (0..dy).flat_map(move |y| (0..dx).map(move |x| UVec3::new(x, y, z))))
}

#[derive(Clone)]
pub struct PlaceMap<'s> {
    pub solvent: &'s Structure,
    dimensions: UVec3,
    // Invariant: length = (solvent.natoms() / 8).ceil() * dimensions.product()
    // Invariant: Any wasted bits are always set to zero.
    // Value is 0 if the place of this solvent bead is not otherwise occupied.
    // Conversely, if the value of a bit is 1, a solvent bead cannot be placed there.
    //
    // To illustrate what these wasted bits are, the figure below shows meaningful bits as `*`, and
    // the wasted bits as `0`. The bytes are stored in first to last order, and the bits in the
    // bytes are just treated according to their significance.
    //
    // ******** ******** 0000****
    //                   ^^^^-- Wasted bits.
    placements: Box<[u8]>,
}

impl<'s> PlaceMap<'s> {
    pub fn new(solvent: &'s Structure, size: UVec3) -> Self {
        let n_cells = size.x as usize * size.y as usize * size.z as usize;
        let n_bytes = solvent.natoms().div_ceil(8);
        let placements = vec![0x00; n_bytes * n_cells].into_boxed_slice();
        Self {
            dimensions: size,
            solvent,
            placements,
        }
    }

    pub fn n_cells(&self) -> usize {
        self.dimensions.x as usize * self.dimensions.y as usize * self.dimensions.z as usize
    }

    pub fn iter_placements(&self) -> impl Iterator<Item = (Vec3, Placement)> {
        iter_3d(self.dimensions).map(|cell_pos| {
            let translation = cell_pos.as_vec3() * self.solvent.boxvecs.as_vec3();
            let placement = self.get(cell_pos).unwrap();
            (translation, placement)
        })
    }

    pub fn iter_placements_mut(&mut self) -> impl Iterator<Item = PlacementMut> {
        let n_beads = self.solvent.natoms();
        let n_bytes = n_beads.div_ceil(8);
        self.placements
            .chunks_exact_mut(n_bytes)
            .map(move |bytes| PlacementMut::new(bytes, n_beads))
    }

    pub fn iter_atoms_chunks(&self) -> impl Iterator<Item = Box<[Atom]>> + '_ {
        self.iter_placements().map(|(translation, placement)| {
            let solvent_positions =
                self.solvent
                    .atoms()
                    .zip(placement)
                    .filter_map(move |(&sb, occupied)| {
                        if occupied {
                            None
                        } else {
                            let mut sb = sb; // Copy the solvent bead.
                            sb.position += translation.convert();
                            Some(sb)
                        }
                    });
            solvent_positions.collect::<Box<[_]>>()
        })
    }

    pub fn iter_atoms(&self) -> impl Iterator<Item = Atom> + '_ {
        self.iter_placements().flat_map(|(translation, placement)| {
            self.solvent
                .atoms()
                .zip(placement)
                .filter_map(move |(&sb, occupied)| {
                    if occupied {
                        None
                    } else {
                        let mut sb = sb; // Copy the solvent bead.
                        sb.position += translation.convert();
                        Some(sb)
                    }
                })
        })
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = eightyseven::structure::Vec3> + '_ {
        self.iter_atoms().map(|atom| atom.position)
    }

    pub fn get_mut(&mut self, pos: UVec3) -> Option<PlacementMut> {
        if pos.x >= self.dimensions.x || pos.y >= self.dimensions.y || pos.z >= self.dimensions.z {
            return None;
        }

        // TODO: Is it wasteful to 'compute' this here? Or does it reduce down well?
        let n_beads = self.solvent.natoms();
        let n_bytes = n_beads.div_ceil(8);
        let idx = index_3d(pos, self.dimensions);
        // This unwrap should be safe, since we checked at the start of the function.
        let bytes = self.placements.chunks_exact_mut(n_bytes).nth(idx).unwrap();

        Some(PlacementMut::new(bytes, n_beads))
    }

    pub fn get(&self, pos: UVec3) -> Option<Placement> {
        if pos.x >= self.dimensions.x || pos.y >= self.dimensions.y || pos.z >= self.dimensions.z {
            return None;
        }

        // TODO: Is it wasteful to 'compute' this here? Or does it reduce down well?
        let n_beads = self.solvent.natoms();
        let n_bytes = n_beads.div_ceil(8);
        let idx = index_3d(pos, self.dimensions);
        // This unwrap should be safe, since we checked at the start of the function.
        let bytes = self.placements.chunks_exact(n_bytes).nth(idx).unwrap();

        Some(Placement::new(bytes, n_beads))
    }

    /// Return the count of the positions that are marked as occupied.
    pub fn occupied_count(&self) -> u64 {
        // Since a bit is marked as 1 iff the position it represents is occupied, we can just
        // perform a popcount over the whole bit array.
        // Note that we don't have to worry about the wasted bits, since these are always set to
        // zero.
        self.placements.iter().map(|b| b.count_ones() as u64).sum()
    }

    /// Return the count of unoccupied positions.
    pub fn unoccupied_count(&self) -> u64 {
        let total = self.solvent.natoms() * self.n_cells();
        total as u64 - self.occupied_count()
    }

    /// Repair all wasted bits in the type.
    fn repair(&mut self) {
        let n_beads = self.solvent.natoms();
        let n_bytes = n_beads.div_ceil(8);
        self.placements
            .chunks_exact_mut(n_bytes)
            .map(|chunk| PlacementMut::new(chunk, n_beads))
            .for_each(|mut p| p.repair())
    }
}

impl std::ops::BitOrAssign for PlaceMap<'_> {
    // Invariant: The or operation will not break the invariant as long as the two operands are
    // also valid.
    fn bitor_assign(&mut self, rhs: Self) {
        assert_eq!(self.dimensions, rhs.dimensions);
        assert_eq!(self.n_cells(), rhs.n_cells());
        assert_eq!(self.solvent, rhs.solvent);

        self.placements
            .iter_mut()
            .zip(rhs.placements)
            .for_each(|(s, r)| {
                *s |= r;
            });
    }
}

impl std::ops::BitAndAssign for PlaceMap<'_> {
    // Invariant: The and operation will not break the invariant as long as the two operands are
    // also valid.
    fn bitand_assign(&mut self, rhs: Self) {
        assert_eq!(self.dimensions, rhs.dimensions);
        assert_eq!(self.n_cells(), rhs.n_cells());
        assert_eq!(self.solvent, rhs.solvent);

        self.placements
            .iter_mut()
            .zip(rhs.placements)
            .for_each(|(s, r)| {
                *s &= r;
            });
    }
}

impl std::ops::BitXorAssign for PlaceMap<'_> {
    // Invariant: The xor operation will not break the invariant as long as the two operands are
    // also valid.
    fn bitxor_assign(&mut self, rhs: Self) {
        assert_eq!(self.dimensions, rhs.dimensions);
        assert_eq!(self.n_cells(), rhs.n_cells());
        assert_eq!(self.solvent, rhs.solvent);

        self.placements
            .iter_mut()
            .zip(rhs.placements)
            .for_each(|(s, r)| {
                *s ^= r;
            });
    }
}

impl std::ops::Not for PlaceMap<'_> {
    type Output = Self;

    fn not(mut self) -> Self::Output {
        // In order to preserve the invariant that the wasted bits are always zero, we need to be
        // careful, here.
        // First, we will go over all the placement bytes and apply a not.
        self.placements.iter_mut().for_each(|s| *s = !(*s));

        // After that, we go and repair the wasted bits, setting them to zero.
        self.repair();

        self
    }
}

pub struct Placement<'ps> {
    pub(crate) bytes: &'ps [u8],
    pub(crate) len: usize,
    idx: usize,
}

impl<'ps> Placement<'ps> {
    fn new(bytes: &'ps [u8], len: usize) -> Self {
        Self { bytes, len, idx: 0 }
    }
}

impl Iterator for Placement<'_> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.len {
            return None;
        }

        let byte = self.bytes[self.idx / 8];
        let bit = (byte >> (self.idx % 8)) & 1;
        self.idx += 1;
        Some(bit > 0)
    }
}

pub struct PlacementMut<'ps> {
    bytes: &'ps mut [u8],
    len: usize,
}

impl<'ps> PlacementMut<'ps> {
    pub(crate) fn new(bytes: &'ps mut [u8], len: usize) -> Self {
        Self { bytes, len }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    /// Set the bit value for a solvent bead to represent that its spot is occupied.
    ///
    /// This communicates that a solvent particle cannot be rendered at this spot.
    pub fn set_occupied(&mut self, idx: usize) {
        if idx >= self.len {
            panic!(
                "index out of range (length was {} but index was {idx})",
                self.len
            )
        }

        let byte = &mut self.bytes[idx / 8];
        *byte |= 1 << (idx % 8);
    }

    /// Get the value of the bit at `idx`.
    ///
    /// # Panics
    ///
    /// This function will panic if the provided `idx` exceeds the range represented by this value.
    pub fn get(&self, idx: usize) -> bool {
        if idx >= self.len {
            panic!(
                "index out of range (length was {} but index was {idx})",
                self.len
            )
        }

        let byte = self.bytes[idx / 8];
        byte >> (idx % 8) & 1 > 0
    }

    /// Repair the wasted bits at the end of this [`PlacementMut`] by setting them to zero.
    ///
    /// This operation is idempotent.
    fn repair(&mut self) {
        let n_bytes = self.bytes.len();
        let n_bits = n_bytes * u8::BITS as usize;
        if n_bits == self.len {
            return;
        }

        let Some(last) = self.bytes.last_mut() else {
            return;
        };
        let n_wasted = n_bits - self.len;

        let mask = (1u8 << (u8::BITS as usize - n_wasted)) - 1;
        *last &= mask;
    }
}

#[cfg(test)]
mod tests {
    use eightyseven::{reader::ReadGro, structure::Structure};

    use super::*;

    #[test]
    fn repair_placement_mut() {
        let mut bytes = [0xff, 0xff, 0xff];
        let mut p = PlacementMut::new(&mut bytes, 20);

        assert_eq!(*p.bytes.last().unwrap(), 0b11111111);
        p.repair();
        assert_eq!(*p.bytes.last().unwrap(), 0b00001111);
        p.repair();
        assert_eq!(*p.bytes.last().unwrap(), 0b00001111);
    }

    /// Verify that the repair that occurs within the not operation works.
    ///
    /// If the repair is incorrect, we would expect to see `unoccupied_count > natoms_total` after
    /// the not.
    #[test]
    fn check_not_repair() {
        let mut solvent = Structure::open_gro("structures/water.gro").unwrap();
        // Make sure we have a structure size that does neatly fit in some number of bytes.
        if solvent.natoms() % 8 == 0 {
            solvent.atoms.pop();
        }
        let natoms = solvent.natoms();
        let size = UVec3::new(3, 5, 7);
        let mut placemap = PlaceMap::new(&solvent, size);
        let natoms_total = (natoms * placemap.n_cells()) as u64;
        let n_occupied = 0;

        // Set some atom to be occupied.
        placemap
            .get_mut(UVec3::new(1, 3, 5))
            .unwrap()
            .set_occupied(161);
        // Update natoms_total and n_occupied to reflect this.
        let natoms_total = natoms_total - 1;
        let n_occupied = n_occupied + 1;

        assert_eq!(placemap.unoccupied_count(), natoms_total);
        assert_eq!(placemap.occupied_count(), n_occupied);

        placemap = !placemap;

        assert_eq!(placemap.unoccupied_count(), n_occupied);
        assert_eq!(placemap.occupied_count(), natoms_total);

        placemap = !placemap;

        assert_eq!(placemap.unoccupied_count(), natoms_total);
        assert_eq!(placemap.occupied_count(), n_occupied);
    }
}
