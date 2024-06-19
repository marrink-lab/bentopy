use std::io;

use clap::Parser;
use glam::Mat3;
use rand::{seq::SliceRandom, Rng};

use args::Args;
use config::Configuration;
use placement::{Batch, Placement, PlacementList};
use state::State;

mod args;
mod config;
mod placement;
mod state;

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
        println!("  idx \tname      \tnrots\ttarget\tplaced\ttime (s)\tremark");
        println!("  ----\t----------\t-----\t------\t------\t--------\t------");
        let mut nrots_tot = 0;
        let mut target_tot = 0;
        let mut hits_tot = 0;
        for (i, (name, nrots, target, hits, duration)) in self.entries.iter().enumerate() {
            let ok = if hits == target { " " } else { "<" };
            println!(
                "{i:>4}\t{name:>10}\t{nrots:>5}\t{target:>6}\t{hits:>6}\t{duration:8.2}\t{ok}"
            );
            nrots_tot += nrots;
            target_tot += target;
            hits_tot += hits;
        }
        println!( "    \t          \t{nrots_tot:>5}\t{target_tot:>6}\t{hits_tot:>6}\t{packing_duration:8.2}")
    }
}

fn main() -> io::Result<()> {
    eprintln!("Welcome to new pack. Let's see if I can finish this by Friday.");

    // Preparation.
    let args = Args::parse();
    let config: Configuration = serde_json::from_reader(std::fs::File::open(&args.path)?)?;
    let mut state = State::new(args, config)?;

    // Packing.
    let mut placements = Vec::new();
    let mut summary = Summary::new();
    let packing_start = std::time::Instant::now();

    for segment in &mut state.segments {
        // Prepare the session.
        let start = std::time::Instant::now();
        let mut placement = Placement::new(segment.name.clone(), segment.path.clone());
        state.space.enter_session(&segment.compartments);

        // Get the free locations.
        // TODO: Putting it all into a Vec seems dumb and a waste. Maybe we can select a random
        // subset if the number of locations is very large?
        let shuffle_guess = segment.number;
        let mut locations = state.space.get_free_locations().collect::<Vec<_>>();
        let (mut shuffled, mut locations) =
            locations.partial_shuffle(&mut state.rng, shuffle_guess);
        let mut cursor = 0;

        let mut hits = 0;
        let mut batch_positions = Vec::new();
        'placement: while hits < segment.number {
            // TODO: This can become more efficient for successive placement failures.
            let voxels = match segment.voxels() {
                Some(voxels) => voxels,
                None => {
                    segment.voxelize(state.space.resolution, state.bead_radius);
                    segment.voxels().unwrap()
                }
            };

            let [maxx, maxy, maxz] = {
                let [bx, by, bz] = state.space.dimensions;
                let [sx, sy, sz] = voxels.dimensions();
                [bx - sx, by - sy, bz - sz]
            };

            // Pick a random location.
            let &location = loop {
                let candidate = match shuffled.get(cursor) {
                    // Get it from the shuffled set.
                    Some(location) => location,
                    // Or, shuffle more positions if available.
                    None if !locations.is_empty() => {
                        (shuffled, locations) =
                            locations.partial_shuffle(&mut state.rng, shuffle_guess);
                        cursor = 0;
                        shuffled.first().unwrap()
                    }
                    None => {
                        // We've gone through all the locations.
                        // TODO: This is rather unlikely, but must be dealt with correctly. I think
                        // this may require a nice custom type that does some internal mutabilty
                        // trickery.
                        // For now, we'll just break here.
                        break 'placement;
                    }
                };
                cursor += 1;
                let &[x, y, z] = candidate;
                if x < maxx && y < maxy && z < maxz {
                    break candidate;
                }
            };

            // Try it.
            if state.space.check_collisions(voxels, location) {
                // We found a good spot. Stomp on those stamps!
                state.space.stamp(voxels, location);
                // Transform the location to nm.
                batch_positions.push(location.map(|v| v as f32 * state.space.resolution));

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
        state.space.exit_session();
        let duration = std::time::Instant::now().duration_since(start);

        // Save the batches.
        summary.push(
            segment.name.clone(),
            placement.n_batches(),
            segment.number,
            hits,
            duration.as_secs_f64(),
        );
        placements.push(placement);
    }

    let packing_duration = std::time::Instant::now()
        .duration_since(packing_start)
        .as_secs_f64();
    eprintln!("Packing process took {packing_duration:.3} s.");

    // Output.
    let placement_list_path = format!(
        "{}_placements.json",
        state.output.dir.join(&state.output.title).to_str().unwrap()
    );
    let placement_list_file = std::fs::File::create(&placement_list_path)?;
    let placement_list = PlacementList::new(placements, &state);
    serde_json::to_writer(placement_list_file, &placement_list)?;
    eprintln!("Wrote placement list to {placement_list_path:?}.");

    // Final summary.
    eprintln!("Summary:");
    summary.present(packing_duration);

    Ok(())
}
