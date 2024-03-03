use std::io;
use std::path::PathBuf;

use clap::Parser;

use crate::solvate::solvate;
use crate::structure::Structure;

mod placement;
mod solvate;
mod structure;

/// Solvate.
#[derive(Debug, Parser)]
struct Args {
    /// Structure input path.
    input: PathBuf,
    /// Solvent template path.
    template: PathBuf,
    /// Output path.
    output: PathBuf,
    /// Cutoff radius.
    #[arg(long, default_value = "0.4")]
    cutoff: f32,
    /// Center the structure in the new box.
    #[arg(short, long)]
    center: bool,
    // /// Number of positive ions.
    // #[arg(short, long, default_value = "0")]
    // posions: u32,
    // /// Number of negative ions.
    // #[arg(short, long, default_value = "0")]
    // negions: u32,

    // TODO: Consider adding an option for overshooting, cutting off, or undershooting when the
    // input box size is not an integer multiple of the template box.
}

fn main() -> io::Result<()> {
    let config = Args::parse();

    let template = Structure::read_from_gro_file(config.template)?;
    let mut structure = Structure::read_from_gro_file(config.input)?;

    solvate(
        &mut structure,
        &template,
        config.cutoff,
        config.center,
        // config.posions,
        // config.negions,
    );

    structure.write_to_gro_file(config.output)
}
