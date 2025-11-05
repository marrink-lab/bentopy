use std::path::PathBuf;

use clap::{Parser, ValueEnum};

/// Provide additional version information, including the current git hash.
mod version {
    pub(crate) struct Version {
        pkg_version: &'static str,
        git_version: &'static str,
    }

    impl Into<clap::builder::Str> for Version {
        fn into(self) -> clap::builder::Str {
            let Self {
                pkg_version,
                git_version,
            } = self;
            format!("{pkg_version} ({git_version})").into()
        }
    }

    pub(crate) const VERSION: Version = Version {
        pkg_version: env!("CARGO_PKG_VERSION"),
        git_version: git_version::git_version!(
            args = ["--broken", "--always", "--exclude", "*"],
            prefix = "git:",
            fallback = "?"
        ),
    };
}

/// Pack a space.
#[derive(Debug, Parser)]
#[command(version = version::VERSION)]
pub struct Args {
    // Configuration input file to define the run (json).
    pub config: PathBuf,

    // Placement list output file (json).
    pub output: PathBuf,

    /// Sort the input structures by approximate size to optimize packing.
    #[arg(long, default_value = "moment-of-inertia")]
    pub rearrange: RearrangeMethod,

    /// Random number generator seed.
    #[arg(long)]
    pub seed: Option<u64>,

    /// Set the bead radius that is considered during voxelization in nm.
    #[arg(long, default_value_t = 0.20)]
    pub bead_radius: f32,

    /// Disable printing the summary after the placement procedure.
    #[arg(long)]
    pub no_summary: bool,

    /// Set the multiplier for determining the maximum number of tries for placing a segment.
    ///
    /// The number of segments requested is multiplied with this value to set the upper limit on
    /// how many unsuccessful attempts will be made to place it.
    ///
    /// Setting a higher value will result in more attempts, which can help improve the total
    /// number of segments packed in very dense systems.
    #[arg(long, default_value_t = 1000)]
    pub max_tries_mult: u64,

    /// Set the divisor for determining the maximum number of tries for placing a segment per
    /// rotation.
    ///
    /// This is a rather advanced option, and it is unlikely that it requires tweaking.
    ///
    /// The maximum number of attempts (see `--max-tries-mult`) is divided by this divisor value to
    /// determine the number of unsuccessful attempts that are allowed for a given rotation. After
    /// this number of unsuccessful tries, the next random rotation will be tried.
    #[arg(long, default_value_t = 100)]
    pub max_tries_rot_div: u64,

    /// Display verbose output.
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum RearrangeMethod {
    MomentOfInertia,
    Volume,
    BoundingSphere,
    /// Keep the arrangement specified in the input file.
    None,
}
