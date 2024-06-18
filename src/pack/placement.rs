use std::path::PathBuf;

use serde::Serialize;

use crate::config::TopolIncludes;
use crate::state::{Position, Size, State};

type Rotation = [[f32; 3]; 3];

#[derive(Serialize)]
pub struct Batch {
    rotation: Rotation,
    positions: Vec<Position>,
}

impl Batch {
    pub fn new(rotation: Rotation, positions: Vec<Position>) -> Self {
        Self {
            rotation,
            positions,
        }
    }
}

#[derive(Serialize)]
pub struct Placement {
    name: String,
    path: PathBuf,
    batches: Vec<Batch>,
}

impl Placement {
    pub fn new(name: String, path: PathBuf) -> Self {
        Self {
            name,
            path,
            batches: Default::default(),
        }
    }

    pub fn push(&mut self, batch: Batch) {
        self.batches.push(batch)
    }

    pub fn n_batches(&self) -> usize {
        self.batches.len()
    }
}

// TODO: These types are prime targets for moving into some `core` crate so that both render and
// this can use the type.
#[derive(Serialize)]
pub struct PlacementList {
    title: String,
    size: Size,
    topol_includes: TopolIncludes,
    placements: Vec<Placement>,
}

impl PlacementList {
    pub fn new(placements: impl IntoIterator<Item = Placement>, state: &State) -> Self {
        Self {
            title: state.output.title.to_string(),
            size: state.space.size,
            topol_includes: state.output.topol_includes.clone().unwrap_or_default(),
            placements: placements.into_iter().collect(),
        }
    }
}
