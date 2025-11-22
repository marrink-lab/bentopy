//! Render a placement list to a gro file with speed.
//!
//! By Marieke Westendorp, 2024.
//! <ma3ke.cyber@gmail.com>
use std::path::PathBuf;

use clap::{Parser, command};

use crate::limits::Limits;
use crate::render::{Mode, ResnumMode, render};

mod limits;
mod render;
mod structure;

/// Render structures from a placement list into a gro file.
///
///
/// Structures specified in the placement list are retrieved from their pdb or gro
/// files and placed into a gro file according to their rotations and positions.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Args {
    /// Path to the placement list.
    ///
    /// To read from stdin, pass "-".
    input: PathBuf,

    /// Output gro file path.
    output: PathBuf,

    /// Write a topology (.top) file.
    #[arg(short, long)]
    topol: Option<PathBuf>,

    /// Root path for the structure paths.
    ///
    /// When set, this path will be prepended to any relative path pointing to a structure in the
    /// placement list. Absolute paths are respected.
    #[arg(long)]
    root: Option<PathBuf>,

    /// Granularity of the produced output.
    #[arg(long, value_enum, default_value_t = Mode::Full, conflicts_with="topol")]
    mode: Mode,

    /// Write out a unique resnum for each segment instance, or use one grouped resnum for each
    /// instance of a segment.
    #[arg(long, value_enum, default_value_t)]
    resnum_mode: ResnumMode,

    /// Only render structures that have a position within a smaller cuboid.
    ///
    /// Arguments can be provided as a comma-separated array of 6 values. Each value can be a
    /// number indicating a bound or a non-numerical value indicating an unset bound.
    #[arg(long)]
    limits: Option<Limits>,

    #[arg(long)]
    ignore_tags: bool,
}

fn main() -> anyhow::Result<()> {
    let Args {
        input,
        output,
        topol,
        root,
        mode,
        limits,
        resnum_mode,
        ignore_tags,
    } = Args::parse();
    render(
        input,
        output,
        topol,
        root,
        limits,
        mode,
        resnum_mode,
        ignore_tags,
    )
}
