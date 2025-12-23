use std::path::PathBuf;
use std::str::FromStr;

use clap::{Parser, ValueEnum};
use eightyseven::structure::ResName;

/// Solvate.
#[derive(Debug, Parser)]
#[command(about, version = bentopy::core::version::VERSION)]
pub struct Args {
    /// Structure input path.
    #[arg(short, long)]
    pub input: PathBuf,
    /// Output path.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Lowest allowable distance between solvent and structure beads (nm).
    ///
    /// This is the minimum allowable center-to-center distance when checking for a collision
    /// between a solvent bead and a structure bead.
    ///
    /// For the solvent-solvent cutoff distance, see `--solvent-cutoff`.
    ///
    /// For Martini solvation, the default cutoff is 0.43 nm. For atomistic solvation, the default
    /// cutoff is 0.28.
    #[arg(long)]
    pub cutoff: Option<f32>,

    /// Lowest allowable distance between solvent beads (nm).
    ///
    /// This is the minimum allowable center-to-center distance between solvent beads when cutting
    /// the sides to fit in the specified output structure box.
    ///
    /// For Martini solvation, the default cutoff is 0.21 nm. For atomistic solvation, the default
    /// cutoff is 0.21.
    #[arg(long)]
    pub solvent_cutoff: Option<f32>,

    /// List of resnames to ignore when checking against structure-solvent collisions.
    #[arg(long, value_delimiter = ',')]
    pub ignore: Vec<ResName>,

    /// The type of water written to the output file.
    #[arg(long, value_enum, default_value_t)]
    pub water_type: WaterType,

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

    /// Set the charge to neutralize with additional ions.
    ///
    /// In essence, this is a shorthand for explicitly providing ions to compensate the charge as
    /// substitutes. This can be helpful for automation purposes.
    ///
    /// By default, NA (positive) and CL (negative) ions are used, but a different
    /// positive-negative substitution pair can be specified by prepending a colon followed by the
    /// positive and negative substitute name separated by one comma.
    ///
    ///     <charge>:<positive ion name>,<negative ion name>
    ///
    /// So, by default, some `<charge>` is interpreted as
    ///
    ///     <charge>:NA,CL
    #[arg(long, allow_hyphen_values = true)]
    pub charge: Option<ChargeConfig>,

    /// Combine substitutes with identical names into one block.
    #[arg(long, default_value_t)]
    pub no_combine_substitutes: bool,

    /// Set whether and in what way substitutes are sorted.
    #[arg(long, value_enum, default_value_t)]
    pub sort_substitutes: SortBehavior,

    /// Random number generator seed for solvent substitution, such as ion placement.
    #[arg(long)]
    pub seed: Option<u64>,

    /// If the solvent template contains velocity information, write these velocities to the output
    /// file.
    #[arg(long, default_value_t)]
    pub write_velocities: bool,

    /// Append solvation topology lines to a path.
    #[arg(short = 't', long)]
    pub append_topol: Option<PathBuf>,

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

#[derive(Debug, Default, Clone, Copy, ValueEnum, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub enum WaterType {
    #[default]
    Martini,
    Tip3P,
}

impl WaterType {
    pub const fn default_cutoff(&self) -> f32 {
        match self {
            WaterType::Martini => 0.43,
            WaterType::Tip3P => 0.28,
        }
    }

    pub const fn default_solvent_cutoff(&self) -> f32 {
        match self {
            WaterType::Martini => 0.21,
            WaterType::Tip3P => 0.21,
        }
    }
}

#[derive(Debug, Default, Clone, ValueEnum, PartialEq, Eq)]
pub enum BoundaryMode {
    /// Cut at the structure box size and remove solvent residues that overlap with the periodic
    /// neighbors.
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
pub struct ChargeConfig {
    charge: i64,
    positive: String,
    negative: String,
}

impl FromStr for ChargeConfig {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split(':');
        let charge = words
            .next()
            .ok_or("expected a charge")?
            .parse::<i64>()
            .map_err(|err| err.to_string())?;
        let (positive, negative) = if let Some(names) = words.next() {
            names
                .split_once(',')
                .ok_or("expected two substitute names separated by a comma")?
        } else {
            ("NA", "CL")
        };
        Ok(Self {
            charge,
            positive: positive.to_string(),
            negative: negative.to_string(),
        })
    }
}

impl ChargeConfig {
    /// Bake the [`ChargeConfig`] into a [`Substitute`] if applicable.
    ///
    /// If the charge is 0, no `Substitute` is returned since there is no charge to neutralize.
    pub fn bake(self) -> Option<Substitute> {
        // Choose the name of the ion to neutralize with. If the charge to compensate is negative,
        // we compensate with the positive ion substitute, and vice versa.
        let name = match self.charge {
            ..0 => self.positive,
            0 => return None,
            1.. => self.negative,
        };

        Some(Substitute {
            name,
            number: self.charge.unsigned_abs(),
        })
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

#[derive(ValueEnum, Default, Clone, Debug)]
pub enum SortBehavior {
    #[default]
    Size,
    RevSize,
    Alphabetical,
    RevAlphabetical,
    No,
}
