use std::{
    io::{self, Read},
    path::Path,
};

use eightyseven::structure::{AtomName, AtomNum, ResName, ResNum};
use glam::{Mat3, Vec3};

/// A molecule as derived from PDB 'ATOM' records.
#[derive(Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl Molecule {
    /// Create a new molecule from a PDB file.
    fn from_pdb(pdb: &str) -> Self {
        let mut atoms = Vec::new();
        for line in pdb.lines() {
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let atom = Atom::from_pdb_atom_line(line);
                atoms.push(atom);
            }
        }

        Self { atoms }
    }

    /// Create a new molecule from a gro file.
    fn from_gro(gro: &str) -> Self {
        let mut atoms = Vec::new();
        let mut lines = gro.lines();
        let _title = lines.next().expect("expected a title line");
        let n_atoms: usize = lines
            .next()
            .expect("expected the number of atoms")
            .trim()
            .parse()
            .expect("could not parse the number of atoms");
        for line in lines.take(n_atoms) {
            let atom = Atom::from_gro_atom_line(line);
            atoms.push(atom);
        }

        // let _boxvec = lines.next().expect("expected the box vectors");

        Self { atoms }
    }

    /// Apply a `tag` to the atoms in this [`Molecule`].
    ///
    /// This will replace the `resname` for each atom with the provided `tag`.
    pub fn apply_tag(&mut self, tag: impl Into<ResName>) {
        let tag = tag.into();
        for Atom { resname, .. } in &mut self.atoms {
            *resname = tag;
        }
    }

    /// Returns the center of this [`Molecule`].
    ///
    /// Here, center refers to mean position. This can be understood as a center of mass where all
    /// atoms are considered to have the same weight.
    pub fn center(&self) -> Vec3 {
        let mut mean = Vec3::ZERO;
        for atom in &self.atoms {
            mean += atom.pos
        }
        mean / self.atoms.len() as f32
    }

    /// Translate all atoms such that their center lies at the origin.
    ///
    /// Here, center refers to mean position. This can be understood as a center of mass where all
    /// atoms are considered to have the same weight.
    pub fn translate_to_center(&mut self) {
        let center = self.center();
        for atom in &mut self.atoms {
            atom.pos -= center;
        }
    }

    /// Returns the minimum _x_, _y_, and _z_ values of this [`Molecule`] as a [`Vec3`].
    ///
    /// Note that this is a virtual point, not an actual point in the atoms point cloud.
    /// In other words, it is not the lowest point in a sorted list of existing positions.
    pub fn min(&self) -> Vec3 {
        let mut min = Vec3::ZERO;
        for atom in &self.atoms {
            let pos = atom.pos;
            if pos.x < min.x {
                min.x = pos.x
            }
            if pos.y < min.y {
                min.y = pos.y
            }
            if pos.z < min.z {
                min.z = pos.z
            }
        }
        min
    }
}

/// Assumes the [Molecule]'s atoms have been centered.
pub fn rotate_molecule(molecule: &Molecule, rotmat: &Mat3) -> Molecule {
    let mut molecule = molecule.clone();
    for atom in &mut molecule.atoms {
        atom.pos = rotmat.mul_vec3(atom.pos)
    }

    molecule
}

/// The representation of an atom.
#[derive(Clone)]
pub struct Atom {
    /// Atom name.
    pub name: AtomName,
    /// Residue name.
    pub resname: ResName,
    /// Residue sequence number.
    pub resnum: ResNum,
    /// Number of the atom within its structure.
    ///
    /// Also understood as its 'serial' number.
    pub num: AtomNum,
    /// Position in nanometers.
    pub pos: Vec3,
}

impl Atom {
    // COLUMNS        DATA  TYPE    FIELD        DEFINITION
    // -------------------------------------------------------------------------------------
    //  1 -  6        Record name   "ATOM  "
    //  7 - 11        Integer       serial       Atom  serial number.
    // 13 - 16        Atom          name         Atom name.
    // 17             Character     altLoc       Alternate location indicator.
    // 18 - 20        Residue name  resName      Residue name.
    // 22             Character     chainID      Chain identifier.
    // 23 - 26        Integer       resSeq       Residue sequence number.
    // 27             AChar         iCode        Code for insertion of residues.
    // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    // 55 - 60        Real(6.2)     occupancy    Occupancy.
    // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    // 77 - 78        LString(2)    element      Element symbol, right-justified.
    // 79 - 80        LString(2)    charge       Charge  on the atom.
    //
    // Examples:
    // ATOM      1  N   ALA A   1      11.104   6.134  -6.504
    // ATOM      1  N   MET     1      48.048  69.220  58.803
    /// Read a single "ATOM" or "HETATM" record from a PDB and return an Atom.
    fn from_pdb_atom_line(line: &str) -> Atom {
        let serial = line[6..11].trim().parse().unwrap();
        let name = line[12..16].trim().into();
        let resname = line[17..20].trim().into();
        let resnum = line[22..26].trim().parse().unwrap();
        let x = line[30..38].trim().parse().unwrap();
        let y = line[38..46].trim().parse().unwrap();
        let z = line[46..54].trim().parse().unwrap();
        Atom {
            name,
            resname,
            resnum,
            num: serial,
            pos: Vec3::new(x, y, z) / 10.0, // Convert from Å to nm.
        }
    }

    // `    2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244`
    /// Read a single atom line from a gro file and return an Atom.
    fn from_gro_atom_line(line: &str) -> Atom {
        let resnum = line[0..5].trim().parse().unwrap();
        let resname = line[5..10].trim().into();
        let name = line[10..15].trim().into(); // Atom name.
        let num = line[15..20].trim().parse().unwrap(); // Atom number.
        let x = line[20..28].trim().parse().unwrap();
        let y = line[28..36].trim().parse().unwrap();
        let z = line[36..44].trim().parse().unwrap();
        Atom {
            name,
            resname,
            resnum,
            num,
            pos: Vec3::new(x, y, z), // Values are already in nm.
        }
    }

    /// Create a dummy atom.
    pub fn dummy() -> Self {
        Self {
            name: "DUMMY".into(),
            resname: "DUMMY".into(),
            resnum: 0,
            num: 0,
            pos: Vec3::ZERO,
        }
    }
}

/// Load a [`Molecule`] from a pdb file.
pub fn load_molecule<P: AsRef<Path> + std::fmt::Debug>(path: P) -> io::Result<Molecule> {
    let mut data = String::new();
    eprintln!("\tLoading {path:?}...");
    std::fs::File::open(&path)?.read_to_string(&mut data)?;

    let molecule = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gro") => Molecule::from_gro(&data),
        Some("pdb") => Molecule::from_pdb(&data),
        None | Some(_) => {
            eprintln!("WARNING: Assuming {path:?} is a pdb file.");
            Molecule::from_pdb(&data)
        }
    };

    Ok(molecule)
}
