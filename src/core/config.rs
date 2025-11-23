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
struct Compartment {
    id: String,
    mask: Mask,
}

#[derive(Debug)]
enum Mask {
    Voxels(PathBuf),
    Shape(Shape),
}

#[derive(Debug)]
// TODO: This name is not appropriate.
enum Shape {
    Sphere { center: Center, radius: f64 },
    Cuboid { start: [f32; 3], end: [f32; 3] },
    Combination { expression: Expression },
}

#[derive(Debug)]
enum Center {
    Center,
    Point(Point),
}

#[derive(Debug)]
struct Segment {
    name: String,
    tag: Option<String>, // Should be 5-character ArrayString.
    path: PathBuf,
    compartment_ids: Box<[String]>,
    quantity: Quantity,
}

#[derive(Debug)]
enum Quantity {
    /// Copy number.
    Number(u64),
    /// Molarity (mol/L).
    Concentration(f64),
}

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
