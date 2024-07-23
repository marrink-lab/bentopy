use std::path::PathBuf;

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

    /// Cutoff radius (nm).
    #[arg(long, default_value_t = 0.21)]
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

    // /// Number of positive ions.
    // #[arg(short, long, default_value = "0")]
    // posions: u32,
    // /// Number of negative ions.
    // #[arg(short, long, default_value = "0")]
    // negions: u32,
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
