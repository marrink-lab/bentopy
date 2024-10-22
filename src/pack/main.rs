use std::io;

use args::Args;
use clap::Parser;
use config::Configuration;
use placement::PlacementList;
use state::State;

mod args;
mod config;
mod mask;
mod placement;
mod rules;
mod session;
mod state;
mod structure;

const CLEAR_LINE: &str = "\u{1b}[2K\r";

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
        println!( "    \t          \t{perc:>5.1}%\t{target_tot:>6}\t{hits_tot:>6}\t{packing_duration:8.2}\t{ok}")
    }
}

fn main() -> io::Result<()> {
    // Preparation.
    let args = Args::parse();
    let config: Configuration = serde_json::from_reader(std::fs::File::open(&args.config)?)
        .unwrap_or_else(|err| {
            eprintln!(
                "ERROR: Failed to process configuration from {:?}: {err}",
                &args.config
            );
            std::process::exit(1)
        });
    let mut state = State::new(args, config).unwrap_or_else(|err| {
        eprintln!("ERROR: Failed to set up program state: {err}");
        std::process::exit(1)
    });

    // Check whether the rules make any sense.
    state.check_rules().unwrap_or_else(|err| {
        eprintln!("ERROR: Encountered a problem while checking the rules: {err}");
        std::process::exit(1)
    });

    // Packing.
    let packing_start = std::time::Instant::now();
    let (placements, summary) = state.pack(&mut io::stderr())?;
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
    let placement_list_file = std::fs::File::create(&placement_list_path)?;
    let placement_list = PlacementList::new(placements, &state);
    serde_json::to_writer(placement_list_file, &placement_list)?;
    eprintln!("Done in {:.3} s.", start_output.elapsed().as_secs_f64());

    Ok(())
}
