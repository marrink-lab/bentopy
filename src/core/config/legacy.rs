use std::{path::PathBuf, str::FromStr};

use serde::Deserialize;

pub use super::compartment_combinations::Expression as CombinationExpression;
use crate::core::config::{Axes, CompartmentID, Dimensions};

impl<'de> Deserialize<'de> for CombinationExpression {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        CombinationExpression::from_str(&s).map_err(serde::de::Error::custom)
    }
}

/// Avogadro's number (per mol).
const N_A: f64 = 6.0221415e23;

// TODO: I think it's cursed that we store these defaults here. I'd like to create a file
// collecting all of these consts, one day.
fn bead_radius_default() -> f32 {
    0.20 // nm
}

fn max_tries_mult_default() -> u64 {
    1000
}

fn max_tries_rot_div_default() -> u64 {
    100
}

#[derive(Deserialize)]
pub struct General {
    pub seed: Option<u64>,
    #[serde(default = "bead_radius_default")]
    pub bead_radius: f32,
    #[serde(default = "max_tries_mult_default")]
    pub max_tries_mult: u64,
    #[serde(default = "max_tries_rot_div_default")]
    pub max_tries_rot_div: u64,
}

impl Default for General {
    fn default() -> Self {
        Self {
            seed: Default::default(),
            bead_radius: bead_radius_default(),
            max_tries_mult: max_tries_mult_default(),
            max_tries_rot_div: max_tries_rot_div_default(),
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Shape {
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
pub enum Mask {
    Shape(Shape),
    Analytical {
        shape: Shape,
        center: Option<[f32; 3]>,
        radius: Option<f32>,
    },
    Voxels {
        path: PathBuf,
    },
    Combination(CombinationExpression),
}

#[derive(Deserialize)]
pub struct Compartment {
    pub id: CompartmentID,
    #[serde(flatten)]
    pub mask: Mask,
}

impl Compartment {
    pub fn is_predefined(&self) -> bool {
        match &self.mask {
            Mask::Shape(_) | Mask::Analytical { .. } | Mask::Voxels { .. } => true,
            Mask::Combination(_) => false,
        }
    }
}

pub(crate) fn true_by_default() -> bool {
    true
}

#[derive(Deserialize)]
pub struct Space {
    pub size: Dimensions,
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
pub struct Config {
    #[serde(default)]
    pub general: General,
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,
}

mod convert {
    use super::*;
    use crate::core::config;

    impl Into<config::Mask> for Mask {
        fn into(self) -> config::Mask {
            match self {
                Mask::Shape(shape) => config::Mask::Shape(match shape {
                    Shape::Spherical => config::Shape::space_filling_sphere(),
                    Shape::Cuboid | Shape::None => config::Shape::space_filling_cuboid(),
                }),
                Mask::Analytical {
                    shape: Shape::Spherical,
                    center,
                    radius,
                } => config::Mask::Shape(config::Shape::Sphere {
                    center: match center {
                        None => config::Center::Center,
                        Some(center) => config::Center::Point(center),
                    },
                    radius: radius.expect("TODO") as f64,
                }),
                Mask::Analytical {
                    shape: Shape::Cuboid | Shape::None,
                    center: _, // TODO: Error if not None.
                    radius: _, // TODO: Error if not None.
                } => config::Mask::Shape(config::Shape::space_filling_cuboid()),
                Mask::Voxels { path } => config::Mask::Voxels(path),
                Mask::Combination(expression) => {
                    // TODO: This conversion sucks and will be changed.
                    config::Mask::Shape(config::Shape::Combination { expression })
                }
            }
        }
    }

    impl Into<config::Compartment> for Compartment {
        fn into(self) -> config::Compartment {
            let Self { id, mask } = self;
            config::Compartment {
                id,
                mask: mask.into(),
            }
        }
    }

    impl Into<config::Segment> for Segment {
        fn into(self) -> config::Segment {
            let Self {
                name,
                tag,
                quantity,
                path,
                compartments,
                rules,
                rotation_axes,
                initial_rotation,
            } = self;
            // TODO: Implement parsing for these segment properties.
            if !rules.is_empty() {
                todo!("implement segment rules for bent")
            }
            if rotation_axes == Default::default() {
                todo!("implement segment rotation axes for bent")
            }
            if initial_rotation == <[f32; 3]>::default() {
                todo!("implement segment initial rotation for bent")
            }
            config::Segment {
                name,
                tag,
                path,
                compartment_ids: compartments.into_boxed_slice(),
                quantity: quantity.into(),
            }
        }
    }

    impl Into<config::Quantity> for Quantity {
        fn into(self) -> config::Quantity {
            match self {
                Quantity::Number(n) => config::Quantity::Number(n as u64),
                Quantity::Concentration(c) => config::Quantity::Concentration(c),
            }
        }
    }

    impl Into<config::Config> for Config {
        fn into(self) -> config::Config {
            let Self {
                general:
                    General {
                        seed,
                        bead_radius,
                        max_tries_mult,
                        max_tries_rot_div,
                    },
                space:
                    Space {
                        size: dimensions,
                        resolution,
                        compartments,
                        periodic,
                    },
                segments,
                output:
                    Output {
                        title,
                        topol_includes,
                    },
            } = self;

            config::Config {
                title,
                seed,
                bead_radius,
                max_tries_mult,
                max_tries_rot_div,
                dimensions,
                resolution,
                periodic,
                compartments: compartments.into_iter().map(Into::into).collect(),
                topol_includes: topol_includes.unwrap_or_default(),
                segments: segments.into_iter().map(Into::into).collect(),
            }
        }
    }
}
