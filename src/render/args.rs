use std::{path::PathBuf, str::FromStr};

use crate::limits::Limits;
use clap::{Parser, ValueEnum, command};

/// Render structures from a placement list into a gro file.
///
///
/// Structures specified in the placement list are retrieved from their pdb or gro
/// files and placed into a gro file according to their rotations and positions.
#[derive(Debug, Parser)]
#[command(about, version = bentopy::core::version::VERSION)]
pub struct Args {
    /// Path to the placement list.
    ///
    /// To read from stdin, pass "-".
    pub input: PathBuf,

    /// Output gro file path.
    pub output: PathBuf,

    /// Write a topology (.top) file.
    #[arg(short, long)]
    pub topol: Option<PathBuf>,

    /// Root path for the structure paths.
    ///
    /// When set, this path will be prepended to any relative path pointing to a structure in the
    /// placement list. Absolute paths are respected.
    #[arg(long)]
    pub root: Option<PathBuf>,

    /// Granularity of the produced output.
    #[arg(long, value_enum, default_value_t = Mode::Full, conflicts_with="topol")]
    pub mode: Mode,

    /// Write out a unique resnum for each segment instance, or use one grouped resnum for each
    /// instance of a segment.
    #[arg(long, value_enum, default_value_t)]
    pub resnum_mode: ResnumMode,

    /// Only render structures that have a position within a smaller cuboid.
    ///
    /// Arguments can be provided as a comma-separated array of 6 values. Each value can be a
    /// number indicating a bound or a non-numerical value indicating an unset bound.
    #[arg(long)]
    pub limits: Option<Limits>,

    #[arg(long)]
    pub ignore_tags: bool,
}

#[derive(Debug, Default, Clone, ValueEnum)]
pub enum Mode {
    /// Output all atoms.
    #[default]
    Full,
    /// Only the backbone atoms.
    Backbone,
    /// Only the alpha carbons.
    Alpha,
    /// One bead per residue.
    Residue,
    /// One bead per instance of a structure.
    Instance,
}

impl FromStr for Mode {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "full" => Self::Full,
            "backbone" | "bb" => Self::Full,
            "alpha" | "ca" | "a" => Self::Alpha,
            "residue" | "res" => Self::Residue,
            "instance" => Self::Instance,
            _ => return Err(()),
        })
    }
}

#[derive(Debug, Default, Clone, ValueEnum)]
pub enum ResnumMode {
    /// Each segment instance has a unique residue number.
    #[default]
    Instance,
    /// All instances of a segment have the same residue number.
    Segment,
}

impl FromStr for ResnumMode {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "instance" => Self::Instance,
            "segment" => Self::Segment,
            _ => return Err(()),
        })
    }
}
