use eightyseven::structure::Atom;
use glam::Vec3;

use crate::convert::Convert;

/// A cuboid water box.
pub struct WaterBox<const NATOMS: usize, const NATOMS_PER_RES: usize> {
    // TODO: Make this ResName and AtomName.
    resname: &'static str,
    names: [&'static str; NATOMS_PER_RES],
    /// Box dimensions.
    dimensions: Vec3,
    /// A list of the positions.
    ///
    /// For atomistic ordering, the canonical water ordering is oxygen first, followed by
    /// hydrogens. That is, the first atom in the residue is the atom that is taken as the position
    /// for the solvent-structure collision checks.
    positions: [Vec3; NATOMS],
}

impl<const NATOMS: usize, const NATOMS_PER_RES: usize> WaterBox<NATOMS, NATOMS_PER_RES> {
    /// Creates a new [`WaterBox`].
    pub const fn new(
        resname: &'static str,
        names: [&'static str; NATOMS_PER_RES],
        dimensions: Vec3,
        positions: [Vec3; NATOMS],
    ) -> Self {
        // The number of atoms must be a multiple of the number of atoms per residue.
        const { assert!(NATOMS.is_multiple_of(NATOMS_PER_RES)) }

        Self {
            resname,
            names,
            dimensions,
            positions,
        }
    }

    /// Returns solvent residue name for this [`WaterBox`].
    pub const fn resname(&self) -> &'static str {
        self.resname
    }

    /// Returns the number of points per residue for this [`WaterBox`].
    pub const fn residue_points(&self) -> usize {
        NATOMS_PER_RES
    }

    /// Returns the dimensions of this [`WaterBox`].
    pub const fn dimensions(&self) -> Vec3 {
        self.dimensions
    }

    /// Return an iterator over the positions must be considered for solvent-structure collisions.
    pub fn positions(&'_ self) -> PositionsIterator<'_> {
        PositionsIterator {
            step: NATOMS_PER_RES,
            idx: 0,
            positions: &self.positions,
        }
    }

    pub fn residues(&'_ self) -> ResiduesIterator<'_> {
        ResiduesIterator {
            step: NATOMS_PER_RES,
            idx: 0,
            resname: self.resname,
            names: &self.names,
            positions: &self.positions,
        }
    }
}

pub struct PositionsIterator<'wb> {
    step: usize,
    idx: usize,
    positions: &'wb [Vec3],
}

impl<'wb> Iterator for PositionsIterator<'wb> {
    type Item = Vec3;

    fn next(&mut self) -> Option<Self::Item> {
        let idx = self.idx;
        self.idx += self.step;
        self.positions.get(idx).copied()
    }
}

pub struct ResiduesIterator<'wb> {
    step: usize,
    idx: usize,
    resname: &'static str,
    names: &'wb [&'static str],
    positions: &'wb [Vec3],
}

pub struct Residue<'wb> {
    resname: &'static str,
    names: &'wb [&'static str],
    positions: &'wb [Vec3],
}

impl<'wb> Residue<'wb> {
    // TODO: Maybe this better fits where it is used. Here we start to introduce implementation
    // details that are more local at the call site.
    pub fn atoms<'a: 'wb>(
        &'wb self,
        resnum: &'a mut u32,
        atomnum: &'a mut u32,
        translation: Vec3,
    ) -> impl Iterator<Item = Atom> + use<'wb> {
        let rn = *resnum;
        *resnum += 1;

        Iterator::zip(self.positions.iter(), self.names).map(move |(&pos, &name)| {
            let an = *atomnum;
            *atomnum += 1;

            Atom {
                resnum: rn,
                resname: self.resname.into(),
                atomname: name.into(),
                atomnum: an,
                position: (pos + translation).convert(),
                velocity: Default::default(),
            }
        })
    }
}

impl<'wb> Iterator for ResiduesIterator<'wb> {
    type Item = Residue<'wb>;

    fn next(&mut self) -> Option<Self::Item> {
        let idx = self.idx;
        self.idx += self.step;
        let positions = self.positions.get(idx..self.idx)?;
        Some(Self::Item {
            resname: self.resname,
            names: self.names,
            positions,
        })
    }
}

// Include the actual water coordinates.
include!("martini.rs");
include!("tip3p.rs");
