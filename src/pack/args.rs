use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use bentopy::core::{citation::CITATION, version::VERSION};

/// Pack a space.
#[derive(Debug, Parser)]
#[command(about, version = VERSION, after_help = CITATION)]
pub struct Args {
    // Configuration input file to define the run (bent).
    pub config: PathBuf,

    // Placement list output file (json).
    pub output: PathBuf,

    /// Sort the input structures by approximate size to optimize packing.
    ///
    /// Different heuristics are available.
    #[arg(long, default_value = "moment")]
    pub rearrange: RearrangeMethod,

    /// Random number generator seed.
    #[arg(long)]
    pub seed: Option<u64>,

    /// Set the bead radius that is considered during voxelization in nm.
    ///
    /// If set, this value overwrites any value provided in the input file.
    /// If this value is not provided in the input file, the default value is 0.20 nm.
    #[arg(long)]
    pub bead_radius: Option<f64>,

    /// Set the multiplier for determining the maximum number of tries for placing a segment.
    ///
    /// The number of segments requested is multiplied with this value to set the upper limit on
    /// how many unsuccessful attempts will be made to place it.
    ///
    /// Setting a higher value will result in more attempts, which can help improve the total
    /// number of segments packed in very dense systems.
    ///
    /// If set, this value overwrites any value provided in the input file.
    /// If this value is not provided in the input file, the default value is 1000.
    #[arg(long)]
    pub max_tries_mult: Option<u64>,

    /// Set the divisor for determining the maximum number of tries for placing a segment per
    /// rotation.
    ///
    /// This is a rather advanced option, and it is unlikely that it requires tweaking.
    ///
    /// The maximum number of attempts (see `--max-tries-mult`) is divided by this divisor value to
    /// determine the number of unsuccessful attempts that are allowed for a given rotation. After
    /// this number of unsuccessful tries, the next random rotation will be tried.
    ///
    /// If set, this value overwrites any value provided in the input file.
    /// If this value is not provided in the input file, the default value is 100.
    #[arg(long)]
    pub max_tries_rot_div: Option<u64>,

    /// Disable printing the summary after the placement procedure.
    #[arg(long)]
    pub no_summary: bool,

    /// Display verbose output.
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum RearrangeMethod {
    /// Use the geometric moment (sum of squared distances from geometric center).
    #[value(alias = "moment-of-inertia")]
    Moment,
    Volume,
    BoundingSphere,
    /// Keep the arrangement specified in the input file.
    None,
}
