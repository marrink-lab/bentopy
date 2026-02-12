use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use eightyseven::structure::ResName;

use bentopy::core::{citation::CITATION, version::VERSION};

use crate::topology::determine_system_charge;
use crate::water::Water;

/// Solvate.
#[derive(Debug, Parser)]
#[command(about, version = VERSION, after_help = CITATION)]
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
    ///
    /// A substitute is defined according to the scheme <name>:<quantity>. For
    /// example, 150 mM NaCl can be described as `-s NA:0.15Ms -s CL:0.15Ms`.
    ///
    /// A shorthand for this is `-s NA,CL:0.15Ms`. Note that this shorthand
    /// respects stoichiometry for ions such as MgCl₂: `-s MG,CL@2:0.026Ms`. The
    /// `@2` here is used to indicate the stoichiometric ratio of Mg:Cl::1:2 for
    /// the dissociated magnesium chloride.
    ///
    /// Quantities can be specified as follows.
    ///
    /// - Molar concentration with respect to solvent quantity:
    ///   floating point number followed by an 'Ms' suffix.
    ///   (Example: `-s NA:0.15Ms` replaces 150 mM of solvent residues with NA,
    ///   determined based on the remaining solvent quantity.)
    ///
    /// - Molar concentration with respect to box volume:
    ///   floating point number followed by an 'M' suffix.
    ///   (Example: `-s NA:0.15M` replaces 150 mM of solvent residues with NA,
    ///   determined based on the box volume.)
    ///   This behaviour is equivalent to that of many other solvation tools, such as `gmx genion`.
    // See https://gitlab.com/gromacs/gromacs/-/blob/release-2026/src/gromacs/gmxpreprocess/genion.cpp#L553
    ///
    /// - Count: an unsigned integer.
    ///   (Example: `-s NA:100` replaces 100 solvent residues with NA.)
    ///
    /// - Ratio: a floating point number.
    ///   (Example: `-s NA:0.1` replaces 10% of solvent residues are replaced with NA.)
    #[arg(short, long = "substitute")]
    pub substitutes: Vec<SubstituteConfig>,

    /// Neutralize the system charge with additional ions.
    ///
    /// By passing `--charge neutral`, the total system charge will be determined automically. This
    /// requires a topology.
    ///
    /// The charge to be neutralized can also be set explicitly by providing an integer.
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
    ///
    /// Here, `<charge>` can be the string 'neutral' or an integer.
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

type SubstituteMultiplier = std::num::NonZeroU64;

#[derive(Debug, Clone)]
pub struct SubstituteConfig {
    names: Box<[(String, SubstituteMultiplier)]>,
    quantity: Quantity,
}

impl SubstituteConfig {
    /// Bake into a [`Substitute`] according to a final volume in nm³ and some number of valid
    /// solvent beads.
    pub fn bake(
        self,
        volume: f64,
        solvent_beads: u64,
        water: Water,
    ) -> impl Iterator<Item = Substitute> {
        let base_number = self.quantity.number(volume, solvent_beads, water);
        self.names.into_iter().map(move |(name, multiplier)| Substitute {
            name,
            number: base_number * multiplier.get(),
        })
    }
}

impl FromStr for SubstituteConfig {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split(':');
        let names = words
            .next()
            .ok_or("expected names, then a colon, followed by a quantity".to_string())?
            .to_string();

        let names = names
            .split(',')
            .map(|name| {
                if let Some((name, ratio)) = name.split_once('@') {
                    let name = name.trim().to_owned();
                    let ratio = parse_ratio(ratio);
                    ratio.map(|ratio: SubstituteMultiplier| (name, ratio))
                } else {
                    Ok((name.trim().to_owned(), 1.try_into().unwrap()))
                }
            })
            .collect::<Result<_, _>>()?;
        let quantity =
            words.next().ok_or("expected a quantifier after the colon".to_string())?.parse()?;
        Ok(Self { names, quantity })
    }
}

fn parse_ratio(ratio: &str) -> Result<SubstituteMultiplier, <SubstituteConfig as FromStr>::Err> {
    ratio
        .parse()
        .map_err(|e| format!("could not parse stoichiometric ratio {ratio:?}: {e}"))
        .and_then(|ratio: u64| {
            SubstituteMultiplier::try_from(ratio).map_err(|e| {
                format!("stoichiometric ratio must be a nonzero unsigned integer: {e}")
            })
        })
}

#[derive(Debug, Clone)]
pub enum Charge {
    Known(i64),
    Neutral,
}

#[derive(Debug, Clone)]
pub struct ChargeConfig {
    charge: Charge,
    positive: String,
    negative: String,
}

impl FromStr for ChargeConfig {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split(':');
        let charge = words.next().ok_or("expected a charge")?;
        let charge = match charge {
            "neutral" => Charge::Neutral,
            n => Charge::Known(n.parse::<i64>().map_err(|err| err.to_string())?),
        };
        let (positive, negative) = if let Some(names) = words.next() {
            names.split_once(',').ok_or("expected two substitute names separated by a comma")?
        } else {
            ("NA", "CL")
        };
        Ok(Self { charge, positive: positive.to_string(), negative: negative.to_string() })
    }
}

impl ChargeConfig {
    /// Bake the [`ChargeConfig`] into a [`Substitute`] if applicable.
    ///
    /// If the charge is 0, no `Substitute` is returned since there is no charge to neutralize.
    pub fn bake<P: AsRef<std::path::Path> + std::fmt::Debug>(
        self,
        topol_path: Option<P>,
    ) -> Result<Option<Substitute>> {
        // Choose the name of the ion to neutralize with. If the charge to compensate is negative,
        // we compensate with the positive ion substitute, and vice versa.
        let charge = match self.charge {
            Charge::Known(known) => known,
            Charge::Neutral => {
                let topol_path =
                    topol_path.context("Automatic charge neutralization requires a topology")?;
                let start = std::time::Instant::now();
                let charge = determine_system_charge(&topol_path).with_context(|| {
                    format!(
                        "Could not determine system charge from topology file at {topol_path:?}"
                    )
                })?;
                let duration = start.elapsed().as_secs_f32();
                eprintln!("According to the topology, the total system charge is {charge}.");
                eprintln!("Determining charges took {duration:.3} s.");
                charge
            }
        };
        let name = match charge {
            ..0 => self.positive,
            0 => return Ok(None),
            1.. => self.negative,
        };

        Ok(Some(Substitute { name, number: charge.unsigned_abs() }))
    }
}

#[derive(Debug, Clone)]
enum Quantity {
    /// A fixed number.
    Number(u64),
    /// A fraction of the valid solvent beads.
    Ratio(f64),
    /// Molarity in mol/L with respect to the box volume.
    MolarityBox(f64),
    /// Molarity in mol/L with respect to the remaining solvent volume.
    MolaritySolvent(f64),
}

impl FromStr for Quantity {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Molarity (with respect to box volume).
        if let Some(molarity) = s.strip_suffix('M') {
            let molarity = molarity.parse::<f64>().map_err(|err| err.to_string())?;
            return Ok(Self::MolarityBox(molarity));
        }

        // Molarity (with respect to remaining solvent volume).
        if let Some(molarity) = s.strip_suffix("Ms") {
            let molarity = molarity.parse::<f64>().map_err(|err| err.to_string())?;
            return Ok(Self::MolaritySolvent(molarity));
        }

        if let Ok(number) = s.parse::<u64>() {
            return Ok(Self::Number(number));
        }

        if let Ok(ratio) = s.parse::<f64>() {
            if !(0.0..=1.0).contains(&ratio) {
                return Err(format!("ratio value must satisfy `1.0 >= r >= 0.0`, found {ratio}"));
            }
            return Ok(Self::Ratio(ratio));
        }

        Err(format!("could not parse quantity {s:?}"))
    }
}

impl Quantity {
    /// Given some volume in nm³ and some number of initial valid solvent beads, determine the
    /// number of beads according to this [`Quantity`].
    fn number(self, volume: f64, solvent_beads: u64, water: Water) -> u64 {
        const N_AVOGADRO: f64 = 6.0221415e23; // per mol.
        let number = match self {
            Quantity::Number(number) => number,
            Quantity::Ratio(ratio) => (solvent_beads as f64 * ratio).round() as u64,
            Quantity::MolarityBox(molarity) => {
                let volume_liter = volume * 1e-24;
                (volume_liter * molarity * N_AVOGADRO).round() as u64
            }
            Quantity::MolaritySolvent(molarity) => {
                /// The molarity of water in mol/L at 1 bar, 300K.
                const MOLARITY_WATER: f64 = 55.34;
                let solvent_quantity = solvent_beads * water.waters_per_res() as u64;
                let n = (molarity * solvent_quantity as f64) / (MOLARITY_WATER + molarity);
                n.round() as u64
            }
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
