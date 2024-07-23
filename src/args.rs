use std::path::PathBuf;
use std::str::FromStr;

use clap::{Parser, ValueEnum};

/// Solvate.
#[derive(Debug, Parser)]
pub struct Args {
    /// Structure input path.
    pub input: PathBuf,
    /// Solvent template path.
    pub template: PathBuf,
    /// Output path.
    pub output: PathBuf,

    /// Cutoff distance (nm).
    ///
    /// This is the minimum allowable center-to-center distance when checking for a collision
    /// between a solvent bead and a structure bead.
    #[arg(long, default_value_t = 0.43)]
    pub cutoff: f32,

    /// Center the structure in the new box.
    ///
    /// Note that this option only has an effect if the boundary mode is set to `grow`.
    #[arg(short, long)]
    pub center: bool,

    /// Set how the boundaries of the box will be treated.
    #[arg(short, long, value_enum, default_value_t)]
    pub boundary_mode: BoundaryMode,

    /// Set how periodic images of the input structure are considered.
    ///
    /// Note that only orthorhombic (cuboid) periodicity is considered, currently.
    // TODO: Consider whether supporting other forms of periodicity is worthwhile.
    #[arg(short, long, value_enum, default_value_t)]
    pub periodic_mode: PeriodicMode,

    /// Define substitutes for solvent residues, such as ions.
    #[arg(short, long = "substitute")]
    pub substitutes: Vec<SubstituteConfig>,

    /// Random number generator seed for solvent substitution, such as ion placement.
    #[arg(long)]
    pub seed: Option<u64>,

    #[arg(long)]
    pub no_write_parallel: bool,
    /// The suggested number of atoms to format at once.
    ///
    /// Setting this to a larger value will allow more atoms to be formatted at once, at the cost
    /// of higher memory consumption. Setting this to a smaller value will lower the memory
    /// footprint at a possible cost of efficiency.
    ///
    /// Note that depending on the number of residues in the template box, the provided value may
    /// be overshot by at most that number.
    #[arg(long, default_value_t = 10000000)]
    pub buffer_size: usize,
}

#[derive(Debug, Default, Clone, ValueEnum, PartialEq, Eq)]
pub enum BoundaryMode {
    /// Cut at the structure box size and remove solvent residues that overlap with the periodic
    /// neigbors.
    #[default]
    Cut,
    /// If necessary, grow the box of the output structure to fit a whole number of template boxes.
    Grow,
}

#[derive(Debug, Default, Clone, ValueEnum, PartialEq, Eq)]
pub enum PeriodicMode {
    /// Include the periodic images of the structure when checking whether a solvent spot is
    /// occupied.
    #[default]
    Periodic,
    /// Ignore the periodic images of the structure for the solvent placement check.
    Ignore,
    /// Treat any structure atoms outside of the output box as an error and exit.
    Deny,
}

#[derive(Debug, Clone)]
pub struct SubstituteConfig {
    name: String,
    quantity: Quantity,
}

impl SubstituteConfig {
    /// Bake into a [`Substitute`] according to a final volume in nm³ and some number of valid
    /// solvent beads.
    pub fn bake(self, volume: f64, solvent_beads: u64) -> Substitute {
        Substitute {
            name: self.name,
            number: self.quantity.number(volume, solvent_beads),
        }
    }
}

impl FromStr for SubstituteConfig {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split(':');
        let name = words
            .next()
            .ok_or("expected a name, then a colon, followed by a quantity".to_string())?
            .to_string();
        let quantity = words
            .next()
            .ok_or("expected a quantifier after the colon".to_string())?
            .parse()?;
        Ok(Self { name, quantity })
    }
}

#[derive(Debug, Clone)]
enum Quantity {
    /// A fixed number.
    Number(u64),
    /// A fraction of the valid solvent beads.
    Ratio(f64),
    /// Molarity in mol/L.
    Molarity(f64),
}

impl FromStr for Quantity {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(molarity) = s.strip_suffix('M') {
            // Molarity.
            return Ok(Self::Molarity(
                molarity.parse::<f64>().map_err(|err| err.to_string())?,
            ));
        }

        if let Ok(number) = s.parse::<u64>() {
            return Ok(Self::Number(number));
        }

        if let Ok(ratio) = s.parse::<f64>() {
            if !(0.0..=1.0).contains(&ratio) {
                return Err(format!(
                    "ratio value must satisfy `1.0 >= r >= 0.0`, found {ratio}"
                ));
            }
            return Ok(Self::Ratio(ratio));
        }

        Err(format!("could not parse quantity {s:?}"))
    }
}

impl Quantity {
    /// Given some volume in nm³ and some number of initial valid solvent beads, determine the
    /// number of beads according to this [`Quantity`].
    fn number(self, volume: f64, solvent_beads: u64) -> u64 {
        const N_AVOGADRO: f64 = 6.0221415e23; // per mol.
        let volume_liter = volume * 1e-24;
        let number = match self {
            Quantity::Number(number) => number,
            Quantity::Ratio(ratio) => (solvent_beads as f64 * ratio).round() as u64,
            Quantity::Molarity(molarity) => (volume_liter * molarity * N_AVOGADRO).round() as u64,
        };
        assert!(
            number <= solvent_beads,
            "the number of substitutions that are specified ({number}) exceeds the number of \
                solvent positions that can be substituted ({solvent_beads})"
        );
        number
    }
}

pub struct Substitute {
    pub name: String,
    pub number: u64,
}
