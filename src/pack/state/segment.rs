use std::path::PathBuf;

use bentopy::core::config::{Axes, CompartmentID, Quantity};
use glam::{Mat3, Quat};

use crate::state::{ORDER, Rotation, Voxels};
use crate::structure::Structure;
use crate::voxelize::voxelize;

pub struct Segment {
    pub name: String,
    pub tag: Option<String>,
    pub quantity: Quantity,
    pub compartments: Box<[CompartmentID]>,
    pub path: PathBuf,
    pub rotation_axes: Axes,
    pub(crate) structure: Structure,
    /// The initial rotation of the structure must be applied before the random rotation.
    pub(crate) initial_rotation: Rotation,
    /// Invariant: This rotation must satisfy the constraints set by the `rotation_axes` field by
    /// construction.
    pub(crate) rotation: Rotation,
    pub(crate) voxels: Option<Voxels>,
}

impl Segment {
    /// Set a new rotation the [`Segment`].
    ///
    /// This invalidates the voxelization.
    ///
    /// The internal `rotation_axes` are taken into account when storing the rotation, such that a
    /// rotation stored in a `Segment` is always internally consistent.
    pub fn set_rotation(&mut self, rotation: Rotation) {
        // FIXME: Assert it's a well-formed rotation?
        // TODO: This seems slightly hacky, since we are converting the rotation between different
        // formats a couple of times. It should be fine---we just lose an unimportant bit of
        // accuracy on a random rotation---but perhaps it is more wise to store the rotation
        // internally as a quaternion and only convert it to Mat3 when writing to the placement
        // list.
        let axes = &self.rotation_axes;
        let (ax, ay, az) = Quat::from_mat3(&rotation).to_euler(ORDER);
        let (ax, ay, az) = (
            if axes.x { ax } else { 0.0 },
            if axes.y { ay } else { 0.0 },
            if axes.z { az } else { 0.0 },
        );
        self.rotation = Mat3::from_euler(ORDER, ax, ay, az);
        self.voxels = None;
    }

    /// Get the correctly formed rotation of this [`Segment`].
    pub fn rotation(&self) -> Rotation {
        self.rotation * self.initial_rotation
    }

    /// Voxelize this [`Segment`] according to its current rotation.
    ///
    /// The voxelization can be accessed through [`Segment::voxels`].
    pub fn voxelize(&mut self, resolution: f32, radius: f32) {
        self.voxels = Some(voxelize(
            &self.structure,
            self.rotation(),
            resolution,
            radius,
        ));
    }

    /// If available, return a reference to the voxels that represent this [`Segment`].
    pub fn voxels(&self) -> Option<&Voxels> {
        self.voxels.as_ref()
    }
}
