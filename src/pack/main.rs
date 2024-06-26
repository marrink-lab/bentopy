use std::io;

use clap::Parser;
use glam::Mat3;
use rand::Rng;

use args::Args;
use config::Configuration;
use placement::{Batch, Placement, PlacementList};
use state::{Locations, State};

mod args;
mod config;
mod mask;
mod placement;
mod state;
mod structure;

const CLEAR_LINE: &str = "\u{1b}[2K\r";

type Location = usize;

struct Summary {
    entries: Vec<(String, usize, usize, usize, f64)>,
}

impl Summary {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    fn push(&mut self, name: String, nrots: usize, target: usize, placed: usize, time: f64) {
        self.entries.push((name, nrots, target, placed, time))
    }

    fn present(&self, packing_duration: f64) {
        println!("idx \tname      \tnrots\ttarget\tplaced\ttime (s)\tremark");
        println!("----\t----------\t-----\t------\t------\t--------\t------");
        let mut nrots_tot = 0;
        let mut target_tot = 0;
        let mut hits_tot = 0;
        for (i, (name, nrots, target, hits, duration)) in self.entries.iter().enumerate() {
            let ok = if hits == target { " " } else { "<" };
            println!(
                "{i:>4}\t{name:<10}\t{nrots:>5}\t{target:>6}\t{hits:>6}\t{duration:8.2}\t{ok}"
            );
            nrots_tot += nrots;
            target_tot += target;
            hits_tot += hits;
        }
        let ok = if hits_tot == target_tot { " " } else { "<" };
        println!( "    \t          \t{nrots_tot:>5}\t{target_tot:>6}\t{hits_tot:>6}\t{packing_duration:8.2}\t{ok}")
    }
}

fn main() -> io::Result<()> {
    // Preparation.
    let args = Args::parse();
    let config: Configuration = serde_json::from_reader(std::fs::File::open(&args.path)?)?;
    let mut state = State::new(args, config)?;
    let mut locations = Locations::new();

    // Packing.
    let mut placements = Vec::new();
    let mut summary = Summary::new();
    let packing_start = std::time::Instant::now();

    let n_segments = state.segments.len();
    for (i, segment) in state.segments.iter_mut().enumerate() {
        if segment.target == 0 {
            eprint!(
                "{prefix}({i:>3}/{n_segments}) Skipping attempt to pack 0 instances of segment '{name}'.{suffix}",
                prefix = if state.verbose { "" } else { CLEAR_LINE },
                i = i + 1,
                name = segment.name,
                suffix = if state.verbose { "\n" } else { "" }
            );
            continue;
        }

        eprint!(
            "{prefix}({i:>3}/{n_segments}) Attempting to pack {target:>5} instances of segment '{name}'.{suffix}",
            prefix = if state.verbose { "" } else { CLEAR_LINE },
            i = i + 1,
            target = segment.target,
            name = segment.name,
            suffix = if state.verbose { "\n" } else { "" }
        );

        // Prepare the session.
        let start_session = std::time::Instant::now();
        let mut session = state.space.enter_session(
            // FIXME: This cloned stuff does not sit right with me.
            segment.compartments.iter().cloned(),
            &mut locations,
            segment.target, // FIXME: Can probably be bigger.
        );

        // Set up the placement record for this segment.
        let mut placement = Placement::new(segment.name.clone(), segment.path.clone());

        let mut hits = 0;
        let mut batch_positions = Vec::new();
        'placement: while hits < segment.target {
            // TODO: This can become more efficient for successive placement failures.
            let voxels = match segment.voxels() {
                Some(voxels) => voxels,
                None => {
                    segment.voxelize(session.resolution(), state.bead_radius);
                    segment.voxels().unwrap()
                }
            };

            let [maxx, maxy, maxz] = {
                let [bx, by, bz] = session.dimensions();
                let [sx, sy, sz] = voxels.dimensions();
                [bx - sx, by - sy, bz - sz]
            };

            // Pick a random location.
            let position = loop {
                let candidate = match session.pop_location(&mut state.rng) {
                    Some(location) => location,
                    None => {
                        // We've gone through all the locations.
                        // TODO: This is rather unlikely, but must be dealt with correctly. I think
                        // this may require a nice custom type that does some internal mutabilty
                        // trickery.
                        // For now, we'll just break here.
                        break 'placement;
                    }
                };
                let position = session.position(candidate).unwrap();
                let [x, y, z] = position;
                if x < maxx && y < maxy && z < maxz {
                    break position;
                }
            };

            // Try it.
            if session.check_collisions(voxels, position) {
                // We found a good spot. Stomp on those stamps!
                session.stamp(voxels, position);
                // Transform the location to nm.
                batch_positions.push(position.map(|v| v as f32 * session.resolution()));

                // Let's write out the batch and rotate the segment again.
                // TODO: Perhaps we'll need a little transpose here.
                let batch = Batch::new(segment.rotation, batch_positions.clone());
                placement.push(batch);
                batch_positions.clear();

                let rotation = Mat3::from_quat(state.rng.gen());
                segment.set_rotation(rotation);

                hits += 1;
            }
        }
        let duration = start_session.elapsed().as_secs_f64();

        if state.verbose {
            let total = packing_start.elapsed().as_secs();
            eprintln!(
                "                      Packed {hits:>5} instances in {duration:6.3} s. [{total} s]",
            );
        }

        // Save the batches.
        summary.push(
            segment.name.clone(),
            placement.n_batches(),
            segment.target,
            hits,
            duration,
        );
        placements.push(placement);
    }

    let packing_duration = packing_start.elapsed().as_secs_f64();
    if !state.verbose {
        eprintln!(); // Go down one line to prevent overwriting the last line.
    }
    eprintln!("Packing process took {packing_duration:.3} s.");

    // Drop this memory hog for good measure.
    drop(locations);

    // Final summary.
    if state.summary {
        summary.present(packing_duration);
    }

    // Output.
    let prefix = state.output.dir.join(&state.output.title);
    let placement_list_path = format!("{}_placements.json", prefix.to_str().unwrap());
    eprint!("Writing placement list to {placement_list_path:?}... ");
    let start_output = std::time::Instant::now();
    let placement_list_file = std::fs::File::create(&placement_list_path)?;
    let placement_list = PlacementList::new(placements, &state);
    serde_json::to_writer(placement_list_file, &placement_list)?;
    eprintln!("Done in {:.3} s.", start_output.elapsed().as_secs_f64());

    Ok(())
}
