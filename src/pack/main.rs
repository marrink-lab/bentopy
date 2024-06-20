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

/// Threshold of spare capacity at which the `locations` [`Vec`] ought to be shrunk.
const SHRINK_THRESHOLD: usize = (100 * 1024_usize.pow(2)) / std::mem::size_of::<state::Position>();

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
    eprintln!("Welcome to new pack. Let's see if I can finish this by Friday.");

    // Preparation.
    let args = Args::parse();
    let config: Configuration = serde_json::from_reader(std::fs::File::open(&args.path)?)?;
    let mut state = State::new(args, config)?;
    let mut locations = Vec::new();

    // Packing.
    let mut placements = Vec::new();
    let mut summary = Summary::new();
    let packing_start = std::time::Instant::now();

    let n_segments = state.segments.len();
    for (i, segment) in state.segments.iter_mut().enumerate() {
        eprintln!(
            "({i:>3}/{n_segments}) Attempting to pack {target:>5} instances of segment '{name}'.",
            target = segment.target,
            name = segment.name
        );

        // Prepare the session.
        let start = std::time::Instant::now();
        // FIXME: This cloned stuff does not sit right with me.
        let mut session = state
            .space
            .enter_session(segment.compartments.iter().cloned());

        // Set up the placement record for this segment.
        let mut placement = Placement::new(segment.name.clone(), segment.path.clone());

        // Get the free locations.
        locations.clear();
        session.get_free_locations(&mut locations);
        if locations.spare_capacity_mut().len() > SHRINK_THRESHOLD {
            // Shrink the `locations` if the spare capacity is getting out of hand.
            locations.shrink_to_fit();
        }
        let shuffle_guess = segment.target;
        let (mut shuffled, mut locations) =
            locations.partial_shuffle(&mut state.rng, shuffle_guess);
        let mut cursor = 0;

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
            if session.check_collisions(voxels, location) {
                // We found a good spot. Stomp on those stamps!
                session.stamp(voxels, location);
                // Transform the location to nm.
                batch_positions.push(location.map(|v| v as f32 * session.resolution()));

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
        let end = std::time::Instant::now();
        let duration = end.duration_since(start).as_secs_f64();

        eprintln!(
            "                      Packed {hits:>5} instances in {duration:6.3} s. [{} s]",
            end.duration_since(packing_start).as_secs()
        );

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
