use std::path::PathBuf;

use clap::Parser;

/// Pack a space.
#[derive(Parser)]
#[command(version)]
pub struct Args {
    // Configuration file to define the run (json).
    pub path: PathBuf,

    /// Sort the input structures by approximate size to optimize packing.
    #[arg(long, default_value_t = true)]
    pub rearrange: bool,

    /// Random number generator seed.
    #[arg(long)]
    pub seed: Option<u64>,

    /// Set the minimum number of random rotations per segment kind.",
    #[arg(long, default_value_t = 10)]
    pub rotations: usize,

    /// Set the bead radius that is considered during voxelization in nm.
    #[arg(long, default_value_t = 0.20)]
    pub bead_radius: f32,

    /// Display verbose output.
    #[arg(short, long)]
    pub verbose: bool,
}
