use std::path::PathBuf;

use glam::Mat3;
use serde::{Deserialize, Serialize};

use crate::core::config::Dimensions;

type Rotation = Mat3;
type Position = [f32; 3];
type RowMajorRotation = [[f32; 3]; 3];

/// A list of instance-based [`Placement`]s of structures in a space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlacementList {
    pub title: String,
    pub size: Dimensions,
    #[serde(flatten)]
    pub meta: Option<Meta>,
    pub topol_includes: Vec<String>,
    pub placements: Vec<Placement>,
}

/// Additional information about the packed structure that is stored in the [`PlacementList`].
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Meta {
    pub seed: u64,
    pub max_tries_mult: u64,
    pub max_tries_per_rotation_divisor: u64,
    pub bead_radius: f32,
}

/// A set of positions associated with some structure.
///
/// The structure is named and is stored at the provided `path`.
///
/// Positions are stored in [`Batch`]es.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Placement {
    pub name: String,
    /// A tag that will replace the associated structure's `resname` field, if present.
    pub tag: Option<String>,
    /// Path to the segment's molecule file.
    pub path: PathBuf,
    pub batches: Vec<Batch>,
}

impl Placement {
    /// Creates a new [`Placement`].
    pub fn new(name: String, tag: Option<String>, path: PathBuf) -> Self {
        Self { name, tag, path, batches: Default::default() }
    }

    /// Returns the optional tag of this [`Placement`].
    pub fn tag(&self) -> Option<&str> {
        self.tag.as_deref()
    }

    /// Push a new [`Batch`] onto this [`Placement`].
    pub fn push(&mut self, batch: Batch) {
        self.batches.push(batch)
    }
}

/// A set of positions that share a specific rotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Batch {
    /// Rotations are stored in row-major order, since they are stored like that in the placement
    /// list.
    pub rotation: RowMajorRotation,
    /// Positions in nm.
    pub positions: Vec<Position>,
}

impl Batch {
    /// Create a new [`Batch`] from a rotation and a set of positions.
    ///
    /// - The locations of the `positions` must be provided in nm. Any resolution adjustments
    ///   must be applied by the caller.
    /// - The provided `rotation` is internally converted and stored in row-major order.
    pub fn new(rotation: Rotation, positions: Vec<Position>) -> Self {
        Self { rotation: rotation.transpose().to_cols_array_2d(), positions }
    }

    /// Returns the 3×3 rotation matrix of this [`Batch`].
    ///
    /// The rotation matrix is returned in column-major order.
    pub fn rotation(&self) -> Mat3 {
        // Because the batch rotation is stored in a row-major order, and the initializer we use
        // here assumes column-major order, we need to transpose the matrix before using it.
        Mat3::from_cols_array_2d(&self.rotation).transpose()
    }
}
