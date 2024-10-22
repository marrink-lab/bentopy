use std::io::{self, Write};

use clap::Parser;
use eightyseven::reader::ReadGro;
use eightyseven::writer::WriteGro;
use rand::SeedableRng;

use crate::args::Args;
pub(crate) use crate::args::{BoundaryMode, PeriodicMode};
use crate::solvate::solvate;
use crate::structure::{write_structure, BoxVecsExtension, Structure};
use crate::substitute::substitute;

mod args;
mod cookies;
mod placement;
mod solvate;
mod structure;
mod substitute;

fn main() -> io::Result<()> {
    let config = Args::parse();

    eprint!("Loading template {:?}... ", config.template);
    let start = std::time::Instant::now();
    let template = Structure::open_gro(config.template)?;
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
    if template.natoms() == 0 {
        eprintln!("ERROR: Template contains no atoms, so solvation cannot proceed.");
        std::process::exit(1);
    }
    eprint!("Loading structure {:?}... ", config.input);
    let start = std::time::Instant::now();
    let mut structure = Structure::open_gro(config.input)?;
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    eprintln!("Solvating...");
    let start = std::time::Instant::now();
    let mut placemap = solvate(
        &mut structure,
        &template,
        config.cutoff,
        config.solvent_cutoff,
        config.center,
        config.boundary_mode,
        config.periodic_mode,
    );
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    let volume = structure
        .boxvecs
        .as_vec3()
        .as_dvec3()
        .to_array()
        .iter()
        .product();
    let substitutes = config
        .substitutes
        .into_iter()
        .map(|sc| sc.bake(volume, placemap.unoccupied_count() as u64))
        .collect::<Vec<_>>();

    let substitutes = if !substitutes.is_empty() {
        eprintln!("Making substitutions...");
        let start = std::time::Instant::now();
        let mut rng = match config.seed {
            Some(seed) => rand::rngs::StdRng::seed_from_u64(seed),
            None => rand::rngs::StdRng::from_entropy(),
        };
        let substitutes = substitute(&mut rng, &mut placemap, &substitutes);
        eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
        substitutes
    } else {
        Vec::new()
    };

    eprintln!("Writing to {:?}...", config.output);
    let start = std::time::Instant::now();
    let file = std::fs::File::create(config.output)?;
    let mut writer = io::BufWriter::new(file);
    let buffer_size = config.buffer_size;
    if config.no_write_parallel {
        write_structure::<false>(
            &mut writer,
            &structure,
            &placemap,
            &substitutes,
            buffer_size,
        )?;
    } else {
        write_structure::<true>(
            &mut writer,
            &structure,
            &placemap,
            &substitutes,
            buffer_size,
        )?;
    }
    eprintln!("Writing took {:.3} s", start.elapsed().as_secs_f32());

    let solvent_name = placemap
        .solvent
        .atoms()
        .next()
        .expect("there is at least one solvent bead")
        .resname
        .as_str();

    // To print out a proper topology, we assume that all solvent beads share the same resname.
    // If these assumptions do not hold, no topology can be printed.
    if !placemap
        .solvent
        .atoms()
        .all(|a| a.resname.as_str() == solvent_name)
    {
        eprintln!(
            "WARNING: Cannot output topology information, since not all beads in the \
            template have the same name (expected {solvent_name:?})."
        );
        std::process::exit(0);
    }

    let mut topology = vec![(solvent_name, placemap.unoccupied_count() as usize)];
    topology.extend(substitutes.iter().map(|s| (s.name(), s.natoms())));
    let mut stdout = io::stdout();
    match &config.append_topol {
        Some(path) => eprintln!("Appending solvent topology lines to {path:?} and stdout:"),
        None => eprintln!("Printing Solvent topology lines to stdout:"),
    }
    let mut topol = config
        .append_topol
        .map(|path| {
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
        })
        .transpose()?;
    for (name, natoms) in topology {
        let s = format!("{name}\t{natoms}\n");
        stdout.write_all(s.as_bytes())?;
        if let Some(topol) = &mut topol {
            topol.write_all(s.as_bytes())?;
        }
    }

    Ok(())
}
