use std::io;

use clap::Parser;
use eightyseven::reader::ReadGro;
use structure::write_structure;

use crate::args::Args;
pub(crate) use crate::args::{BoundaryMode, PeriodicMode};
use crate::solvate::solvate;
use crate::structure::Structure;

mod args;
mod placement;
mod solvate;
mod structure;

fn main() -> io::Result<()> {
    let config = Args::parse();

    eprint!("Loading template {:?}... ", config.template);
    let start = std::time::Instant::now();
    let template = Structure::open_gro(config.template)?;
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
    eprint!("Loading structure {:?}... ", config.input);
    let start = std::time::Instant::now();
    let mut structure = Structure::open_gro(config.input)?;
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    eprintln!("Solvating...");
    let start = std::time::Instant::now();
    let placemap = solvate(
        &mut structure,
        &template,
        config.cutoff,
        config.center,
        config.boundary_mode,
        config.periodic_mode,
    );
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    eprintln!("Writing to {:?}...", config.output);
    let start = std::time::Instant::now();
    let file = std::fs::File::create(config.output)?;
    let mut writer = io::BufWriter::new(file);
    let buffer_size = config.buffer_size;
    if config.no_write_parallel {
        write_structure::<false>(&mut writer, &structure, &template, &placemap, buffer_size)?;
    } else {
        write_structure::<true>(&mut writer, &structure, &template, &placemap, buffer_size)?;
    }
    eprintln!("Writing took {:.3} s", start.elapsed().as_secs_f32());

    Ok(())
}
