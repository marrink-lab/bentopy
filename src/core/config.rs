use anyhow::Result;

use std::path::PathBuf;

use compartment_combinations::Expression;

mod bent;
mod compartment_combinations;
pub mod legacy;

pub type CompartmentID = String;
pub type Dimensions = [f32; 3];
type Point = [f32; 3];

#[derive(Debug)]
pub struct Config {
    // General.
    title: String,
    seed: Option<u64>,
    bead_radius: f32,
    max_tries_mult: u64,
    max_tries_rot_div: u64,
    // TODO: Add RearrangeMethod.
    // Space
    dimensions: Dimensions,
    resolution: f32,
    periodic: bool,
    // Compartments
    compartments: Vec<Compartment>,
    // Includes
    topol_includes: Vec<String>,
    // Segments
    segments: Vec<Segment>,
}

impl Config {
    /// Parse a `.bent` input file.
    pub fn parse_bent(s: &str) -> Result<Self> {
        bent::parse_config(s)
    }

    /// Parse a legacy `.json` input file.
    pub fn parse_legacy_json(_s: &str) -> Result<Self> {
        todo!()
    }
}

#[derive(Debug)]
pub struct Compartment {
    pub id: String,
    pub mask: Mask,
}

#[derive(Debug)]
pub enum Mask {
    Voxels(PathBuf),
    Shape(Shape),
}

#[derive(Debug)]
// TODO: This name is not appropriate.
pub enum Shape {
    Sphere { center: Center, radius: f64 },
    Cuboid { start: [f32; 3], end: [f32; 3] },
    Combination { expression: Expression },
}

impl Shape {
    fn space_filling_cuboid() -> Self {
        Self::Cuboid {
            start: Default::default(),
            end: todo!(),
        }
    }

    fn space_filling_sphere() -> Shape {
        Self::Sphere {
            center: Center::Center,
            radius: todo!(),
        }
    }
}

#[derive(Debug)]
pub enum Center {
    Center,
    Point(Point),
}

#[derive(Debug)]
pub struct Segment {
    pub name: String,
    pub tag: Option<String>, // Should be 5-character ArrayString.
    pub path: PathBuf,
    pub compartment_ids: Box<[String]>,
    pub quantity: Quantity,
}

#[derive(Debug)]
pub enum Quantity {
    /// Copy number.
    Number(u64),
    /// Molarity (mol/L).
    Concentration(f64),
}

#[derive(Debug, PartialEq, Eq)]
pub struct Axes {
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl Default for Axes {
    fn default() -> Self {
        Self {
            x: true,
            y: true,
            z: true,
        }
    }
}

impl std::str::FromStr for Axes {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let len = s.len();
        if len > 3 {
            return Err(format!(
                "an axes string may consist of at most 3 characters, '{s}' has {len} characters"
            ));
        }

        Ok(Self {
            x: s.contains('x'),
            y: s.contains('y'),
            z: s.contains('z'),
        })
    }
}
