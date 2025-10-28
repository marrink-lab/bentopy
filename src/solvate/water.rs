use clap::ValueEnum;
use eightyseven::structure::{Atom, AtomName, AtomNum, ResName, ResNum, Vec3};

#[derive(Debug, Default, Clone, Copy, ValueEnum, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub enum WaterType {
    #[default]
    Martini,
    Tip3P,
}

impl WaterType {
    // TODO: SmallVec optimization?
    pub fn expand(&self, position: Vec3) -> Vec<Atom> {
        fn build_atom(
            resnum: ResNum,
            resname: impl Into<ResName>,
            atomname: impl Into<AtomName>,
            atomnum: AtomNum,
            position: Vec3,
        ) -> Atom {
            Atom {
                resnum,
                resname: resname.into(),
                atomname: atomname.into(),
                atomnum,
                position,
                ..Default::default()
            }
        }

        match self {
            Self::Martini => vec![build_atom(0, "W", "W", 0, position)],
            Self::Tip3P => TIP3P_ATOMS
                .map(|(atomnum, atomname, resname, resnum, pos)| {
                    build_atom(resnum, resname, atomname, atomnum, pos + position)
                })
                .to_vec(),
        }
    }

    pub const fn resname(&self) -> &str {
        match self {
            Self::Martini => "W",
            Self::Tip3P => "SOL",
        }
    }

    pub const fn nresidues(&self) -> usize {
        match self {
            Self::Martini => 1,
            Self::Tip3P => 4,
        }
    }

    const fn residue_points(&self) -> usize {
        match self {
            Self::Martini => 1,
            Self::Tip3P => 3,
        }
    }

    pub const fn total_natoms(&self) -> usize {
        self.residue_points() * self.nresidues()
    }
}

const TIP3P_ATOMS: [(ResNum, &str, &str, AtomNum, Vec3); 12] = [
    (1, "OW", "SOL", 1, Vec3::new(-0.09525, -0.0945, -0.08175)),
    (2, "HW1", "SOL", 1, Vec3::new(-0.12825, -0.0415, -0.00975)),
    (3, "HW2", "SOL", 1, Vec3::new(-0.02225, -0.0435, -0.11775)),
    (4, "OW", "SOL", 2, Vec3::new(-0.08625, 0.0885, 0.09925)),
    (5, "HW1", "SOL", 2, Vec3::new(-0.04725, 0.1145, 0.01525)),
    (6, "HW2", "SOL", 2, Vec3::new(-0.04125, 0.0085, 0.12425)),
    (7, "OW", "SOL", 3, Vec3::new(0.07875, 0.0945, -0.10375)),
    (8, "HW1", "SOL", 3, Vec3::new(0.09675, 0.0335, -0.03275)),
    (9, "HW2", "SOL", 3, Vec3::new(0.16475, 0.1265, -0.13075)),
    (10, "OW", "SOL", 4, Vec3::new(0.10275, -0.0885, 0.08625)),
    (11, "HW1", "SOL", 4, Vec3::new(0.02675, -0.1105, 0.03225)),
    (12, "HW2", "SOL", 4, Vec3::new(0.12275, -0.1695, 0.13425)),
];
