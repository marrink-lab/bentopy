use std::path::PathBuf;

use serde::Serialize;

use crate::config::TopolIncludes;
use crate::state::{Size, State};

type Rotation = [[f32; 3]; 3];
type Position = [usize; 3];

#[derive(Serialize)]
struct Batch {
    rotation: Rotation,
    positions: Vec<Position>,
}

#[derive(Serialize)]
pub struct Placement {
    name: String,
    path: PathBuf,
    batches: Vec<Batch>,
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
