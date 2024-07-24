use std::path::PathBuf;

use clap::{Parser, ValueEnum};

/// Pack a space.
#[derive(Debug, Parser)]
#[command(version)]
pub struct Args {
    // Configuration input file to define the run (json).
    pub config: PathBuf,

    // Placement list output file (json).
    pub output: PathBuf,

    /// Sort the input structures by approximate size to optimize packing.
    #[arg(long, default_value = "moment-of-inertia")]
    pub rearrange: Option<RearrangeMethod>,

    /// Random number generator seed.
    #[arg(long)]
    pub seed: Option<u64>,

    /// Set the bead radius that is considered during voxelization in nm.
    #[arg(long, default_value_t = 0.20)]
    pub bead_radius: f32,

    /// Disable printing the summary after the placement procedure.
    #[arg(long)]
    pub no_summary: bool,

    /// Display verbose output.
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum RearrangeMethod {
    Volume,
    MomentOfInertia,
}
