use std::io;

pub(crate) use glam::Mat3;
use rand::Rng as _;

use crate::placement::{Batch, Placement};
use crate::session::Locations;
use crate::state::State;
use crate::{CLEAR_LINE, Summary};

impl State {
    pub fn pack(&mut self, log: &mut impl io::Write) -> anyhow::Result<(Vec<Placement>, Summary)> {
        let start = std::time::Instant::now();
        let mut locations = Locations::new();
        let mut placements = Vec::new();
        let mut summary = Summary::new();

        let state = self;
        let n_segments = state.segments.len();
        for (i, segment) in state.segments.iter_mut().enumerate() {
            if segment.quantity.is_zero() {
                write!(
                    log,
                    "{prefix}({i:>3}/{n_segments}) Skipping attempt to pack 0 instances of segment '{name}'.{suffix}",
                    prefix = if state.verbose { "" } else { CLEAR_LINE },
                    i = i + 1,
                    name = segment.name,
                    suffix = if state.verbose { "\n" } else { "" }
                )?;
                continue;
            }

            write!(
                log,
                "{prefix}({i:>3}/{n_segments}) Attempting to pack {target:>5} instances of segment '{name}'.{suffix}",
                prefix = if state.verbose { "" } else { CLEAR_LINE },
                i = i + 1,
                target = segment.quantity,
                name = segment.name,
                suffix = if state.verbose { "\n" } else { "" }
            )?;

            // Prepare the session.
            let start_session = std::time::Instant::now();
            let mut session = state.space.enter_session(
                // FIXME: This cloned stuff does not sit right with me.
                segment.compartments.iter().cloned(),
                segment.rules.iter().cloned(),
                &mut locations,
                segment.quantity,
            );

            // Set up the placement record for this segment.
            let mut placement = Placement::new(
                segment.name.clone(),
                segment.tag.clone(),
                segment.path.clone(),
            );

            let mut hits = 0;
            let mut tries = 0; // The number of unsuccessful tries.
            let mut tries_per_rotation = 0;
            // The maximum number of unsuccessful tries.
            let max_tries = state.general.max_tries_multiplier * session.target() as u64;
            // The number of
            let max_tries_per_rotation = max_tries / state.general.max_tries_per_rotation_divisor;
            let mut batch_positions = Vec::new();
            'placement: while hits < session.target() {
                if tries >= max_tries {
                    if state.verbose {
                        writeln!(log, "Exiting after {tries} unsuccessful tries.")?;
                    }
                    break 'placement; // "When you try your best, but you don't succeed."
                }

                // TODO: This can become more efficient for successive placement failures.
                // TODO: Also, this should become a method on Segment.
                let resolution = session.resolution();
                let voxels = match segment.voxels() {
                    Some(voxels) => voxels,
                    None => {
                        segment.voxelize(resolution, state.general.bead_radius);
                        segment.voxels().unwrap() // We just voxelized.
                    }
                };

                // FIXME: Do all this math with glam::U64Vec?
                let [bx, by, bz] = session.dimensions();
                let [sx, sy, sz] = voxels.dimensions();
                if sx > bx || sy > by || sz > bz {
                    tries += 1;

                    // What about another rotation?
                    let rotation = Mat3::from_quat(state.rng.random());
                    segment.set_rotation(rotation);
                    tries_per_rotation = 0; // Reset the counter.
                    if state.verbose {
                        eprintln!("\tNew rotation. The previous rotation would never have fit.")
                    }

                    continue 'placement; // Segment voxels size exceeds the size of the background.
                }
                let [maxx, maxy, maxz] = [bx - sx, by - sy, bz - sz];

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

                    // If we are not packing in a periodic manner, we skip candidate positions that
                    // would necessitate periodic behavior. When a location passes this point, it is
                    // valid for non-periodic placement.
                    if !session.periodic() {
                        // Check if the segment will fit in the background, given this position.
                        let [x, y, z] = position;
                        if x >= maxx || y >= maxy || z >= maxz {
                            continue 'location;
                        }
                    }

                    // All seems good, so we continue to see if this position will be accepted.
                    break position;
                };

                // Reject if it would cause a collision.
                if !session.check_collisions(voxels, position) {
                    tries += 1;
                    tries_per_rotation += 1;
                    if tries_per_rotation >= max_tries_per_rotation {
                        let rotation = Mat3::from_quat(state.rng.random());
                        segment.set_rotation(rotation);
                        tries_per_rotation = 0; // Reset the counter.
                    }
                    continue 'placement; // Reject due to collision.
                }

                // We found a good spot. Stomp on those stamps!
                session.stamp(voxels, position);
                // Transform the location to nm.
                batch_positions.push(position.map(|v| v as f32 * session.resolution()));

                // Let's write out the batch and rotate the segment again.
                // TODO: Perhaps we'll need a little transpose here.
                let batch = Batch::new(segment.rotation(), batch_positions.clone());
                placement.push(batch);
                batch_positions.clear();

                let rotation = Mat3::from_quat(state.rng.random());
                segment.set_rotation(rotation);
                tries_per_rotation = 0; // Reset the counter.

                hits += 1;
            }
            let duration = start_session.elapsed().as_secs_f64();

            if state.verbose {
                let total = start.elapsed().as_secs();
                writeln!(
                    log,
                    "                      Packed {hits:>5} instances in {duration:6.3} s. [{total} s] <{:.4} of max_tries, {:.4} of target>",
                    tries as f32 / max_tries as f32,
                    tries as f32 / session.target() as f32
                )?;
            }

            // Save the batches.
            summary.push(segment.name.clone(), session.target(), hits, duration);
            placements.push(placement);
        }

        if !state.verbose {
            writeln!(log)?; // Go down one line to prevent overwriting the last line.
        }

        Ok((placements, summary))
    }
}
