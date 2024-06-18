use std::path::PathBuf;

use serde::Deserialize;

use crate::state::{CompartmentID, Size};

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum Shape {
    Spherical,
    Cuboid,
    None,
}

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum Mask {
    Shape(Shape),
    Voxels { path: PathBuf },
}

#[derive(Deserialize)]
pub struct Compartment {
    pub id: CompartmentID,
    #[serde(flatten)]
    pub mask: Mask,
}

#[derive(Deserialize)]
pub struct Space {
    pub size: Size,
    pub resolution: f32,
    pub compartments: Vec<Compartment>,
    // TODO: constraint system.
}

#[derive(Deserialize)]
pub struct Segment {
    pub name: String,
    pub number: usize,
    pub path: PathBuf,
    pub compartments: Vec<CompartmentID>,
    // TODO: rotation_axes, center.
}

pub type TopolIncludes = Vec<String>;

#[derive(Deserialize)]
pub struct Output {
    pub title: String,
    pub dir: PathBuf,
    pub topol_includes: Option<TopolIncludes>,
}

#[derive(Deserialize)]
pub struct Configuration {
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,
}
