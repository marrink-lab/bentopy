use clap::ValueEnum;
use eightyseven::structure::Atom;
use glam::Vec3;

use crate::args::WaterType;
use crate::water::models::ResiduesIterator;

mod models;

#[derive(Debug, Default, Clone, Copy, ValueEnum, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub enum Water {
    #[default]
    Martini,
    Tip3P,
}

impl From<WaterType> for Water {
    fn from(value: WaterType) -> Self {
        match value {
            WaterType::Martini => Self::Martini,
            WaterType::Tip3P => Self::Tip3P,
        }
    }
}

impl Water {
    /// Returns an iterator over the positions that must be considered for solvent-structure
    /// collisions.
    pub fn positions(&self) -> impl Iterator<Item = Vec3> {
        match self {
            Water::Martini => models::MARTINI.positions(),
            Water::Tip3P => models::TIP3P.positions(),
        }
    }

    // TODO: All of these delegate functions can be resolved after the refactor.
    /// Returns an iterator over the residues in this [`Water`].
    fn residues(&'_ self) -> ResiduesIterator<'_> {
        match self {
            Water::Martini => models::MARTINI.residues(),
            Water::Tip3P => models::TIP3P.residues(),
        }
    }

    pub const fn dimensions(&self) -> Vec3 {
        match self {
            Water::Martini => models::MARTINI.dimensions(),
            Water::Tip3P => models::TIP3P.dimensions(),
        }
    }

    pub const fn resname(&self) -> &str {
        match self {
            Water::Martini => models::MARTINI.resname(),
            Water::Tip3P => models::TIP3P.resname(),
        }
    }

    // TODO: Rename to natoms_per_residue?
    pub const fn residue_points(&self) -> usize {
        match self {
            Water::Martini => models::MARTINI.residue_points(),
            Water::Tip3P => models::TIP3P.residue_points(),
        }
    }

    /// Returns an iterator over the accepted atoms, ready to be written out.
    pub fn spray<'w, 'a: 'w, A: Iterator<Item = bool>>(
        &'w self,
        mut occupied: A,
        // resnum: &'a mut u32,
        // atomnum: &'a mut u32,
        translation: Vec3,
    ) -> impl Iterator<Item = Atom> + use<'w, A> {
        let mut atomnum = 1;
        let mut resnum = 1;
        // TODO: Add the whole resnum atomnum thing back in!
        self.residues()
            .map(move |res| {
                res.atoms(&mut atomnum, &mut resnum, translation)
                    .collect::<Box<[_]>>()
            })
            .filter(move |_atom| {
                // TODO: Revisit the expect here. Can this happen in normal use? If so, can we make
                // that impossible far before this function is executed. Otherwise, can we make it
                // a nicer error? Or even an unimplemented!().
                !occupied
                    .next()
                    .expect("occupied should have the same length as the number of residues")
            })
            .flatten()
    }

    /// Returns an iterator over the accepted substitute positions.
    pub fn substitute_positions<'w, 'a: 'w, A: Iterator<Item = bool>>(
        &'w self,
        mut occupied: A,
        // resnum: &'a mut u32,
        // atomnum: &'a mut u32,
        translation: Vec3,
    ) -> impl Iterator<Item = Vec3> + use<'w, A> {
        // TODO: Add the whole resnum atomnum thing back in!
        self.positions()
            .map(move |pos| pos + translation)
            .filter(move |_atom| {
                // TODO: Revisit the expect here. Can this happen in normal use? If so, can we make
                // that impossible far before this function is executed. Otherwise, can we make it
                // a nicer error? Or even an unimplemented!().
                !occupied
                    .next()
                    .expect("occupied should have the same length as the number of residues")
            })
    }
}
