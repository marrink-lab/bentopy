use std::path::PathBuf;

use serde::Deserialize;

use crate::state::{Axes, CompartmentID, Size};

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum Shape {
    Spherical,
    Cuboid,
    None,
}

impl std::fmt::Display for Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Shape::Spherical => "spherical",
            Shape::Cuboid => "cuboid",
            Shape::None => "empty ('none')",
        }
        .fmt(f)
    }
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

pub(crate) fn true_by_default() -> bool {
    true
}

#[derive(Deserialize)]
pub struct Space {
    pub size: Size,
    pub resolution: f32,
    pub compartments: Vec<Compartment>,
    #[serde(default = "true_by_default")]
    pub periodic: bool,
    // TODO: constraint system (satisfied _somewhat_ by the notion of a rule).
}

#[derive(Deserialize)]
#[serde(untagged)]
pub enum RuleExpression {
    Rule(String),
    Or(Vec<RuleExpression>),
}

#[derive(Deserialize)]
pub struct Segment {
    pub name: String,
    pub number: usize,
    pub path: PathBuf,
    pub compartments: Vec<CompartmentID>,
    #[serde(default)]
    pub rules: Vec<RuleExpression>,
    #[serde(default)]
    pub rotation_axes: Axes,
    #[serde(default)]
    pub initial_rotation: [f32; 3],
    // TODO: center?
}

pub type TopolIncludes = Vec<String>;

#[derive(Deserialize)]
pub struct Output {
    pub title: String,
    pub topol_includes: Option<TopolIncludes>,
}

#[derive(Deserialize)]
pub struct Configuration {
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,
}
