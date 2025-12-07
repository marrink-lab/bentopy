use std::io::{self, Read};

use anyhow::{Context, bail};
use args::Args;
use ariadne::{Color, Fmt};
use clap::Parser;

use bentopy::core::config::{Config, bent};
use bentopy::core::utilities::CLEAR_LINE;
use state::State;

mod args;
mod mask;
mod session;
mod state;
mod structure;
mod voxelize;

type Location = usize;

struct Summary {
    entries: Vec<(String, usize, usize, f64)>,
}

impl Summary {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    fn push(&mut self, name: String, target: usize, placed: usize, time: f64) {
        self.entries.push((name, target, placed, time))
    }

    fn present(&self, packing_duration: f64) {
        fn percentage(hits: usize, target: usize) -> f32 {
            match target {
                0 => 0.0,
                _ => (hits as f32 / target as f32) * 100.0,
            }
        }

        println!("idx \tname      \t   ok%\ttarget\tplaced\ttime (s)\tremark");
        println!("----\t----------\t------\t------\t------\t--------\t------");
        let mut target_tot = 0;
        let mut hits_tot = 0;
        for (i, (name, target, hits, duration)) in self.entries.iter().enumerate() {
            let perc = percentage(*hits, *target);
            let ok = if hits == target { " " } else { "<" };
            println!(
                "{i:>4}\t{name:<10}\t{perc:>5.1}%\t{target:>6}\t{hits:>6}\t{duration:8.2}\t{ok}"
            );
            target_tot += target;
            hits_tot += hits;
        }
        let perc = percentage(hits_tot, target_tot);
        let ok = if hits_tot == target_tot { " " } else { "<" };
        println!(
            "    \t          \t{perc:>5.1}%\t{target_tot:>6}\t{hits_tot:>6}\t{packing_duration:8.2}\t{ok}"
        )
    }
}

fn main() -> anyhow::Result<()> {
    // Preparation.
    let args = Args::parse();
    let mut state = {
        let config_path = &args.config;
        if config_path.extension().is_some_and(|ext| ext == "json") {
            eprintln!("NOTE: Legacy bentopy input files can be converted to .bent files:");
            let old_path = config_path.to_str().unwrap().fg(Color::Yellow);
            let bent_path = config_path.with_extension("bent");
            let bent_path = bent_path.to_str().unwrap().fg(Color::Cyan);
            eprintln!("You can convert {old_path} to a bent file using the following command:");
            eprintln!();
            eprintln!("\tbentopy-init convert -i {old_path} -o {bent_path}");
            eprintln!();
            bail!("Legacy .json bentopy input files are deprecated")
        }
        let mut file = std::fs::File::open(config_path)
            .with_context(|| format!("Failed to open the configuration file {config_path:?}"))?;
        let mut buf = String::new();
        file.read_to_string(&mut buf)?;
        let config: Config = bent::parse_config(config_path.to_str().unwrap(), &buf)
            .with_context(|| format!("Failed to process configuration from {config_path:?}"))?;

        State::new(args, config).context("Failed to set up program state")?
    };

    // Check whether the compartment masks make any sense.
    state
        .check_masks()
        .with_context(|| format!("Encountered a problem while checking the compartment masks"))?;

    // Packing.
    let packing_start = std::time::Instant::now();
    let (placements, summary) = state
        .pack(&mut io::stderr())
        .context("Encountered a problem while packing")?;
    let packing_duration = packing_start.elapsed().as_secs_f64();
    eprintln!("Packing process took {packing_duration:.3} s.");

    // Final summary.
    if state.summary {
        summary.present(packing_duration);
    }

    // Output.
    let placement_list_path = &state.output.path;
    eprint!("Writing placement list to {placement_list_path:?}... ");
    let start_output = std::time::Instant::now();
    let placement_list_file = std::fs::File::create(&placement_list_path)
        .with_context(|| format!("Failed to create placement list file {placement_list_path:?}"))?;
    let placement_list = state.placement_list(placements);
    serde_json::to_writer(placement_list_file, &placement_list)
        .context("Encountered a problem while writing the placement list")?;
    eprintln!("Done in {:.3} s.", start_output.elapsed().as_secs_f64());

    Ok(())
}
