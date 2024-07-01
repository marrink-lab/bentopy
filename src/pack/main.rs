use std::io;

use clap::Parser;
use glam::Mat3;
use rand::Rng;

use args::Args;
use config::Configuration;
use placement::{Batch, Placement, PlacementList};
use rules::Rule;
use state::{Locations, State};

mod args;
mod config;
mod mask;
mod placement;
#[allow(dead_code)] // FIXME: Remove once the rule system is further developed.
mod rules;
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

    // TODO: Develop rule system further. For now we have an empty rule set.
    let rules = Vec::new();

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
        let mut tries = 0; // The number of unsuccessful tries.
        let max_tries = 1000 * segment.target; // The number of unsuccessful tries.
        let mut batch_positions = Vec::new();
        'placement: while hits < segment.target {
            if tries >= max_tries {
                if state.verbose {
                    eprintln!("Exiting after {tries} unsuccessful tries.");
                }
                break 'placement; // "When you try your best, but you don't succeed."
            }

            let (pre_rules, post_rules) = rules::split(&rules);
            // These rules may be checked many times while picking a valid random location.
            let pre_rules: Vec<&Rule> = pre_rules.collect();

            // TODO: This can become more efficient for successive placement failures.
            // TODO: Also, this should become a method on Segment.
            let resolution = session.resolution();
            let voxels = match segment.voxels() {
                Some(voxels) => voxels,
                None => {
                    segment.voxelize(resolution, state.bead_radius);
                    segment.voxels().unwrap()
                }
            };

            // FIXME: Do all this math with glam::U64Vec?
            let [maxx, maxy, maxz] = {
                let [bx, by, bz] = session.dimensions();
                let [sx, sy, sz] = voxels.dimensions();
                [bx - sx, by - sy, bz - sz]
            };

            // Pick a random location.
            let position = 'location: loop {
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

                // Convert the linear index location to a spatial index.
                let position = session.position(candidate).unwrap();
                let [x, y, z] = position;

                // Check if the segment will fit in the background, given this position.
                if x >= maxx || y >= maxy || z >= maxz {
                    continue 'location;
                }

                // Check the position against the lightweight rules for this segment.
                // FIXME: The math of dimensions considered as the max is kind of shaky I guess.
                for rule in &pre_rules {
                    if !rule.is_satisfied(position, resolution, voxels, session.compartments()) {
                        continue 'location;
                    }
                }

                // All seems good, so we continue to see if this position will be accepted.
                break position;
            };

            // Check if our rules are satisfied.
            for rule in post_rules {
                if !rule.is_satisfied(position, resolution, voxels, session.compartments()) {
                    tries += 1;
                    continue 'placement; // Reject due to breaking a rule.
                }
            }

            // Reject if it would cause a collision.
            if !session.check_collisions(voxels, position) {
                tries += 1;
                continue 'placement; // Reject due to collision.
            }

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
        let duration = start_session.elapsed().as_secs_f64();

        if state.verbose {
            let total = packing_start.elapsed().as_secs();
            eprintln!(
                "                      Packed {hits:>5} instances in {duration:6.3} s. [{total} s] <{:.4} of max_tries, {:.4} of target>",
                tries as f32 / max_tries as f32,
                tries as f32 / segment.target as f32
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
