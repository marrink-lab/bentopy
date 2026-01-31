use std::io::{BufRead, BufReader};
use std::num::NonZeroUsize;

use eightyseven::reader::{ParseList, ReadGro};
use eightyseven::structure::{AtomName, AtomNum, ResName, ResNum};
use eightyseven::writer::{WriteGro, format_atom_line};
use glam::Vec3;

pub type Atom = Vec3;

/// Load a [`Structure`] from a structure file.
///
/// This function will return an error if the structure is empty. Downstream functions may assume
/// that a [`Structure`] has at least one atom.
pub fn load_molecule<P: AsRef<std::path::Path> + std::fmt::Debug>(
    path: P,
) -> anyhow::Result<Structure> {
    use anyhow::{Context, bail};

    let file = std::fs::File::open(&path)?;

    let structure = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gro") => Structure::read_from_file(file)
            .with_context(|| format!("Failed to parse gro file {path:?}"))?,
        Some("pdb") => Structure::read_from_pdb_file(file)
            .with_context(|| format!("Failed to parse PDB file {path:?}"))?,
        None | Some(_) => {
            eprintln!("WARNING: Assuming {path:?} is a PDB file.");
            Structure::read_from_pdb_file(file)
                .with_context(|| format!("Failed to parse the file {path:?} as PDB"))?
        }
    };

    if structure.natoms() == 0 {
        bail!("Structure from {path:?} contains no atoms")
    }

    Ok(structure)
}

/// A structure type that stores its atoms as simple positions.
///
/// Invariant: A `Structure` is always centered, such that its geometric center lies at the origin.
///
/// Invariant: A `Structure` has at least one atom.
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
        // We can safely unwrap because this is the only value we expect.
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

#[non_exhaustive]
#[derive(Debug)]
pub enum ParsePdbError {
    IOError(std::io::Error),
    BadPositionValue {
        err: std::num::ParseFloatError,
        /// Line number.
        ln: NonZeroUsize,
    },
}

impl std::error::Error for ParsePdbError {}

impl std::fmt::Display for ParsePdbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IOError(err) => write!(f, "IO error: {err}"),
            Self::BadPositionValue { err, ln } => {
                write!(f, "Could not parse record on line {ln}: {err}")
            }
        }
    }
}

impl From<std::io::Error> for ParsePdbError {
    fn from(err: std::io::Error) -> Self {
        Self::IOError(err)
    }
}

impl Structure {
    pub fn read_from_pdb_file(file: std::fs::File) -> Result<Structure, ParsePdbError> {
        let reader = BufReader::new(file);

        let mut atoms = Vec::new();
        for (ln, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let ln = NonZeroUsize::new(ln + 1).unwrap(); // We know ln + 1 >= 1 because ln >= 0.
                let bad_position = |err| ParsePdbError::BadPositionValue { err, ln };
                let x = line[30..38].trim().parse().map_err(bad_position)?;
                let y = line[38..46].trim().parse().map_err(bad_position)?;
                let z = line[46..54].trim().parse().map_err(bad_position)?;
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
        // Invariant: A Structure has at least one atom.
        let center = self.atoms().sum::<Atom>() / self.natoms() as f32;
        for pos in &mut self.atoms {
            *pos -= center;
        }
    }

    /// Calculate the moment of inertia for this [`Structure`].
    ///
    /// Invariant: Assumes that the structure is centered.
    ///
    ///     I = Σ(m_i * r_i²)
    pub fn moment_of_inertia(&self) -> f32 {
        self.atoms().map(|atom| atom.length_squared()).sum()
    }

    /// Returns the radius of the point that is farthest from the structure geometric center.
    ///
    /// Invariant: Assumes that the structure is centered.
    pub fn bounding_sphere(&self) -> f32 {
        self.atoms().map(|atom| atom.length()).max_by(f32::total_cmp).unwrap() // Invariant: A structure has at least one atom.
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
