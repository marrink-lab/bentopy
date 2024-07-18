use std::io;
use std::path::PathBuf;

use clap::Parser;
use eightyseven::reader::ReadGro;
use structure::WriteParExt;

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

    /// Cutoff radius (nm).
    #[arg(long, default_value_t = 0.21)]
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

    eprint!("Loading template {:?}... ", config.template);
    let start = std::time::Instant::now();
    let template = Structure::open_gro(config.template)?;
    let delta = std::time::Instant::now() - start;
    eprintln!("Took {:.6} s", delta.as_secs_f32());
    eprint!("Loading structure {:?}... ", config.input);
    let start = std::time::Instant::now();
    let mut structure = Structure::open_gro(config.input)?;
    let delta = std::time::Instant::now() - start;
    eprintln!("Took {:.3} s", delta.as_secs_f32());

    eprint!("Solvating... ");
    let start = std::time::Instant::now();
    solvate(
        &mut structure,
        &template,
        config.cutoff,
        config.center,
        // config.posions,
        // config.negions,
    );
    let delta = std::time::Instant::now() - start;
    eprintln!("Took {:.3} s", delta.as_secs_f32());

    eprintln!("Writing to {:?}...", config.output);
    let start = std::time::Instant::now();
    structure.save_gro_par(config.output)?;
    let delta = std::time::Instant::now() - start;
    eprintln!("Took {:.3} s", delta.as_secs_f32());

    Ok(())
}
