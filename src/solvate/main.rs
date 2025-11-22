use std::cmp::Reverse;
use std::io::{self, Write};

use anyhow::Context;
use args::SortBehavior;
use clap::Parser;
use eightyseven::reader::ReadGro;
use rand::SeedableRng;

use crate::args::Args;
pub(crate) use crate::args::{BoundaryMode, PeriodicMode};
use crate::solvate::solvate;
use crate::structure::{BoxVecsExtension, Structure, write_structure};
use crate::substitute::substitute;
use crate::water::WaterType;

mod args;
mod convert;
mod cookies;
mod placement;
mod solvate;
mod structure;
mod substitute;
mod water;

fn main() -> anyhow::Result<()> {
    let config = Args::parse();

    eprint!("Loading template {:?}... ", config.template);
    let start = std::time::Instant::now();
    let template_path = &config.template;
    let template = {
        let mut template = Structure::open_gro(template_path)
            .with_context(|| format!("Failed to open template structure {template_path:?}"))?;
        if !config.write_velocities {
            // We set all velocities to zero, to indicate that we don't want them written to the
            // output file.
            for atom in &mut template.atoms {
                atom.velocity *= 0.0;
            }
        }
        template
    };
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
    if template.natoms() == 0 {
        eprintln!("ERROR: Template contains no beads, so solvation cannot proceed.");
        std::process::exit(1);
    }

    // We always expect a box of Martini waters.
    let expected_resname = WaterType::Martini.resname();
    if template
        .atoms
        .iter()
        .any(|atom| &atom.resname != expected_resname)
    {
        eprintln!("ERROR: All template beads must be Martini waters ({expected_resname:?}).");
        eprintln!("       The beads in the water template are all treated as the same kind.");
        std::process::exit(1);
    }

    eprint!("Loading structure {:?}... ", config.input);
    let start = std::time::Instant::now();
    let input_path = &config.input;
    let mut structure = Structure::open_gro(input_path)
        .with_context(|| format!("Failed to open structure {input_path:?}"))?;
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    eprintln!("Solvating...");
    let start = std::time::Instant::now();
    let mut placemap = solvate(
        &mut structure,
        &template,
        config.cutoff,
        config.solvent_cutoff,
        &config.ignore,
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
    let neutralizing_ions = config.charge.map(|c| c.bake()).flatten();
    let substitutes = config
        .substitutes
        .into_iter()
        .map(|sc| sc.bake(volume, placemap.unoccupied_count() as u64))
        .chain(neutralizing_ions)
        .collect::<Vec<_>>();

    let substitutions = if !substitutes.is_empty() {
        eprintln!("Making substitutions...");
        let start = std::time::Instant::now();
        let mut rng = match config.seed {
            Some(seed) => rand::rngs::StdRng::seed_from_u64(seed),
            None => rand::rngs::StdRng::from_os_rng(),
        };

        let mut subs = substitute(&mut rng, &mut placemap, &substitutes);

        // If desired, glue the substitutes with identical names together.
        // NOTE: Sorry for the double negative here ;)
        if !config.no_combine_substitutes {
            let mut piles: Vec<substitute::Substitution<'_>> = Vec::new();
            for substitution in subs {
                if let Some(pile) = piles
                    .iter_mut()
                    .find(|pile| pile.name() == substitution.name())
                {
                    // Tack this substitution's replacements onto the existing pile.
                    pile.glue(&substitution)
                } else {
                    piles.push(substitution);
                }
            }

            subs = piles;
        }

        // Order them.
        match config.sort_substitutes {
            SortBehavior::Size => subs.sort_by_cached_key(|s| Reverse(s.natoms())),
            SortBehavior::RevSize => subs.sort_by_cached_key(|s| s.natoms()),
            SortBehavior::Alphabetical => subs.sort_by_key(|s| s.name().to_string()),
            SortBehavior::RevAlphabetical => subs.sort_by_key(|s| Reverse(s.name().to_string())),
            SortBehavior::No => {} // Nothing to do.
        }

        eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
        subs
    } else {
        Vec::new()
    };

    eprintln!("Writing to {:?}...", config.output);
    let start = std::time::Instant::now();
    let output_path = &config.output;
    let file = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create output file {output_path:?}"))?;
    let mut writer = io::BufWriter::new(file);
    let buffer_size = config.buffer_size;
    if config.no_write_parallel {
        write_structure::<false>(
            &mut writer,
            &structure,
            &placemap,
            &substitutions,
            config.water_type,
            buffer_size,
        )
        .with_context(|| {
            format!("Encountered a problem while writing the output structure to {output_path:?}")
        })?;
    } else {
        write_structure::<true>(
            &mut writer,
            &structure,
            &placemap,
            &substitutions,
            config.water_type,
            buffer_size,
        )
        .with_context(|| {
            format!("Encountered a problem while (in parallel) writing the output structure to {output_path:?}")
        })?;
    }
    eprintln!("Writing took {:.3} s", start.elapsed().as_secs_f32());

    let solvent_name = config.water_type.resname();
    let nsolvent = placemap.unoccupied_count() as usize * config.water_type.nresidues();
    let mut topology = vec![(solvent_name, nsolvent)];
    topology.extend(substitutions.iter().map(|s| (s.name(), s.natoms())));
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
                .open(&path)
                .with_context(|| format!("Failed to open topology file {path:?}"))
        })
        .transpose()?;
    for (name, natoms) in topology {
        let s = format!("{name}\t{natoms}\n");
        stdout.write_all(s.as_bytes())?;
        if let Some(topol) = &mut topol {
            topol
                .write_all(s.as_bytes())
                .with_context(|| format!("Failed to append component '{name}' to topology file"))?;
        }
    }

    Ok(())
}
