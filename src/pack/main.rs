use std::io;

use args::Args;
use clap::Parser;
use config::Configuration;
use glam::U64Vec3;
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
        for (i, (name, target, hits, nrots, duration)) in self.entries.iter().enumerate() {
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
        let mut placement = Placement::new(segment.name, segment.path);
        state.space.enter_session(&segment.compartments);

        // Get the free locations.
        let locations = state.space.get_free_locations();

        let mut hits = 0;
        let mut batch_positions = Vec::new();
        while hits < segment.number {
            // TODO: This can become more efficient for successive placement failures.
            let voxels = segment.get_voxels(state.space.resolution, state.bead_radius);
            let u = voxels.indices_where::<true>().map(U64Vec3::from_array);

            // Pick a random location.
            let location = todo!();

            // Try it.
            let l = U64Vec3::from_array(location);
            let v = u.map(|p| p + l);
            for p in v.map(|p| p.to_array()) {}
            if state.space.check_collisions(voxels, location) {
                // We found a good spot. Stomp on those stamps!
                state.space.stamp(voxels, location);
                batch_positions.push(location);

                // Let's write out the batch and rotate the segment again.
                // TODO: Perhaps we'll need a little transpose here.
                let batch = Batch::new(
                    segment.rotation().to_cols_array_2d(),
                    batch_positions.clone(),
                );
                placement.push(batch);
                batch_positions.clear();

                let rotation = todo!("get a random location"); // TODO: This will be a bloody mess again, won't it?
                segment.set_rotation(rotation);

                hits += 1;
            }
        }
        let duration = std::time::Instant::now().duration_since(start);

        // Save the batches.
        placements.push(placement);
        summary.push(
            segment.name,
            placement.n_batches(),
            segment.number,
            hits,
            duration.as_secs_f64(),
        )
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
