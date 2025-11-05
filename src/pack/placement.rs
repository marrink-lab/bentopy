use std::path::PathBuf;

use glam::Mat3;
use serde::Serialize;

use crate::config::TopolIncludes;
use crate::state::{Size, State};

type Rotation = Mat3;
type Position = [f32; 3];
type RowMajorRotation = [[f32; 3]; 3];

#[derive(Serialize)]
pub struct Batch {
    /// Rotations are stored in row-major order, since they are stored like that in the placement
    /// list.
    rotation: RowMajorRotation,
    /// Positions in nm.
    positions: Vec<Position>,
}

impl Batch {
    /// Create a new [`Batch`] from a rotation and a set of positions.
    ///
    /// - The locations of the `positions` must be provided in nm. Any resolution adjustments must
    ///   be applied by the caller.
    /// - The provided `rotation` is internally converted and stored in row-major order.
    pub fn new(rotation: Rotation, positions: Vec<Position>) -> Self {
        Self {
            rotation: rotation.transpose().to_cols_array_2d(),
            positions,
        }
    }
}

#[derive(Serialize)]
pub struct Placement {
    name: String,
    tag: Option<String>,
    path: PathBuf,
    batches: Vec<Batch>,
}

impl Placement {
    pub fn new(name: String, tag: Option<String>, path: PathBuf) -> Self {
        Self {
            name,
            tag,
            path,
            batches: Default::default(),
        }
    }

    pub fn push(&mut self, batch: Batch) {
        self.batches.push(batch)
    }
}

// TODO: These types are prime targets for moving into some `core` crate so that both render and
// this can use the type.
#[derive(Serialize)]
pub struct PlacementList {
    title: String,
    size: Size,
    seed: u64,
    max_tries_mult: u64,
    max_tries_per_rotation_divisor: u64,
    bead_radius: f32,
    topol_includes: TopolIncludes,
    placements: Vec<Placement>,
}

impl PlacementList {
    pub fn new(placements: impl IntoIterator<Item = Placement>, state: &State) -> Self {
        Self {
            title: state.output.title.to_string(),
            size: state.space.size,
            seed: state.seed,
            max_tries_mult: state.max_tries_multiplier,
            max_tries_per_rotation_divisor: state.max_tries_per_rotation_divisor,
            bead_radius: state.bead_radius,
            topol_includes: state.output.topol_includes.clone().unwrap_or_default(),
            placements: placements.into_iter().collect(),
        }
    }
}
