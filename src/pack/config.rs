use std::path::PathBuf;

use serde::Deserialize;

use crate::state::{Axes, CompartmentID, Size};

/// Avogadro's number (per mol).
const N_A: f64 = 6.0221415e23;

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
    Analytical {
        shape: Shape,
        center: Option<[f32; 3]>,
        radius: Option<f32>,
    },
    Voxels {
        path: PathBuf,
    },
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

fn parse_axes<'de, D>(deserializer: D) -> Result<Axes, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.parse().map_err(serde::de::Error::custom)
}

#[derive(Deserialize)]
pub struct Segment {
    pub name: String,
    pub tag: Option<String>,
    #[serde(flatten)]
    pub quantity: Quantity,
    pub path: PathBuf,
    pub compartments: Vec<CompartmentID>,
    #[serde(default)]
    pub rules: Vec<RuleExpression>,
    #[serde(default, deserialize_with = "parse_axes")]
    pub rotation_axes: Axes,
    #[serde(default)]
    pub initial_rotation: [f32; 3],
    // TODO: center?
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "lowercase")]
pub enum Quantity {
    Number(usize),
    /// Concentration in mol/L.
    Concentration(f64),
}

impl Quantity {
    /// Determine the number of segments that is implied by this [`Quantity`].
    ///
    /// In case this `Quantity` is a [`Quantity::Concentration`], the number of segments is
    /// lazily determined from the provided `volume`, and rounded.
    ///
    /// The value returned by `volume` must be in cubic nanometers (nm³).
    pub fn bake<F: Fn() -> f64>(&self, volume: F) -> usize {
        match *self {
            Quantity::Number(n) => n,
            Quantity::Concentration(c) => {
                // n = N_A * c * V
                let v = volume() * 1e-24; // From nm³ to L.
                let n = N_A * c * v;
                f64::round(n) as usize
            }
        }
    }

    /// Returns whether the contained value can be interpreted as resulting in zero placements.
    ///
    /// When the quantity is a `Number(0)` or `Concentration(0.0)`, the baked number is certainly
    /// zero. When `Number(n)` for `n > 0`, the baked number is certainly not zero.
    ///
    /// But, in case of a positive concentration, whether the final number is zero or not depends
    /// on the associated volume.
    ///
    /// If the concentration is smaller than zero, it is treated as a zero.
    pub fn is_zero(&self) -> bool {
        match *self {
            Quantity::Number(n) => n == 0,
            Quantity::Concentration(c) => c <= 0.0,
        }
    }
}

impl std::fmt::Display for Quantity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Quantity::Number(n) => write!(f, "{n} instances"),
            Quantity::Concentration(c) => write!(f, "{c} mol/L"),
        }
    }
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
