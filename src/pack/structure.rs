use std::io::{BufRead, BufReader};

use eightyseven::reader::{ParseList, ReadGro};
use eightyseven::structure::{AtomName, AtomNum, ResName, ResNum};
use eightyseven::writer::{format_atom_line, WriteGro};
use glam::Vec3;

pub type Atom = Vec3;

/// A structure type that stores its atoms as simple positions.
///
/// Invariant: A structure is always centered, such that its geometric center lies at the origin.
pub struct Structure {
    atoms: Vec<Atom>,
}

impl ReadGro<Atom> for Structure {
    const PARSE_LIST: ParseList = ParseList {
        resnum: false,
        resname: false,
        atomname: false,
        atomnum: false,
        position: true,
        velocity: false,
    };

    fn build_atom(
        _resnum: Option<ResNum>,
        _resname: Option<ResName>,
        _atomname: Option<AtomName>,
        _atomnum: Option<AtomNum>,
        position: Option<[f32; 3]>,
        _velocity: Option<[f32; 3]>,
    ) -> Atom {
        Atom::from_array(position.unwrap())
    }

    fn build_structure(
        _title: String,
        atoms: Vec<Atom>,
        _boxvecs: eightyseven::structure::BoxVecs,
    ) -> Self {
        let mut structure = Self { atoms };
        structure.translate_to_center();
        structure
    }
}

impl Structure {
    pub fn read_from_pdb_file(
        file: std::fs::File,
    ) -> Result<Structure, eightyseven::reader::ParseGroError> {
        let reader = BufReader::new(file);

        let mut atoms = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let x = line[30..38].trim().parse().unwrap();
                let y = line[38..46].trim().parse().unwrap();
                let z = line[46..54].trim().parse().unwrap();
                let atom = Atom::new(x, y, z) / 10.0; // Convert from Å to nm.
                atoms.push(atom);
            }
        }

        let mut structure = Self { atoms };
        structure.translate_to_center();
        Ok(structure)
    }

    /// Translate this [`Structure`] such that it's geometric center lies at the origin.
    fn translate_to_center(&mut self) {
        let center = self.atoms().sum::<Atom>() / self.natoms() as f32;
        for pos in &mut self.atoms {
            *pos -= center;
        }
    }

    /// Calculate the moment of inertia for this [`Structure`].
    ///
    /// Assumes that the invariant that the structure has been centered holds.
    ///
    ///     I = Σ(m_i * r_i²)
    pub fn moment_of_inertia(&self) -> f32 {
        self.atoms().map(|atom| atom.length_squared()).sum()
    }
}

impl<'atoms> WriteGro<'atoms, Atom> for Structure {
    fn title(&self) -> String {
        "debug".to_string()
    }

    fn natoms(&self) -> usize {
        self.atoms.len()
    }

    fn atoms(&'atoms self) -> impl Iterator<Item = &'atoms Atom> {
        self.atoms.iter()
    }

    fn boxvecs(&self) -> String {
        "400.0 400.0 400.0".to_string()
    }

    fn format_atom_line(atom: &Atom) -> String {
        format_atom_line(1, "DUMMY", "DUMMY", 2, atom.to_array(), None)
    }
}
