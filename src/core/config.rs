use anyhow::Result;

use std::path::PathBuf;

pub mod bent;
mod compartment_combinations;
pub mod legacy;

pub mod defaults {
    /// Bead radius in nm.
    pub const BEAD_RADIUS: f64 = 0.20;
    pub const MAX_TRIES_MULT: u64 = 1000;
    pub const MAX_TRIES_ROT_DIV: u64 = 100;
}

pub type CompartmentID = String;
pub type Dimensions = [f32; 3];
pub type Point = [f32; 3];

impl Config {
    /// Parse a `.bent` input file.
    pub fn parse_bent(path: &str, s: &str) -> Result<Self> {
        bent::parse_config(path, s)
    }

    /// Parse a legacy `.json` input file.
    pub fn parse_legacy_json(_s: &str) -> Result<Self> {
        todo!()
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub enum RearrangeMethod {
    #[default]
    Moment,
    Volume,
    BoundingSphere,
    None,
}

#[derive(Debug, Default, PartialEq)]
pub struct General {
    pub title: Option<String>,
    pub seed: Option<u64>,
    pub bead_radius: Option<f64>,
    pub max_tries_mult: Option<u64>,
    pub max_tries_rot_div: Option<u64>,
    pub rearrange_method: Option<RearrangeMethod>,
}

#[derive(Debug, Default, PartialEq)]
pub struct Space {
    pub dimensions: Option<Dimensions>,
    pub resolution: Option<f64>,
    pub periodic: Option<bool>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub enum Expr<T> {
    Term(T),
    Not(Box<Self>),
    Or(Box<Self>, Box<Self>),
    And(Box<Self>, Box<Self>),
}

impl<T: std::fmt::Display> std::fmt::Display for Expr<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Term(t) => t.fmt(f),
            Self::Not(expr) => write!(f, "not {expr}"),
            Self::Or(l, r) => write!(f, "({l} or {r})"),
            Self::And(l, r) => write!(f, "({l} and {r})"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl TryFrom<char> for Axis {
    type Error = &'static str;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'x' => Ok(Self::X),
            'y' => Ok(Self::Y),
            'z' => Ok(Self::Z),
            _ => return Err("unknown axis"),
        }
    }
}

impl std::fmt::Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::X => 'x',
            Self::Y => 'y',
            Self::Z => 'z',
        }
        .fmt(f)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    SmallerThan,
    GreaterThan,
}

impl Op {
    fn reverse(self) -> Self {
        match self {
            Self::SmallerThan => Self::GreaterThan,
            Self::GreaterThan => Self::SmallerThan,
        }
    }
}

impl TryFrom<char> for Op {
    type Error = &'static str;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '<' => Ok(Self::SmallerThan),
            '>' => Ok(Self::GreaterThan),
            _ => return Err("unknown operator"),
        }
    }
}

impl std::fmt::Display for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SmallerThan => '<',
            Self::GreaterThan => '>',
        }
        .fmt(f)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Limit {
    pub axis: Axis,
    pub op: Op,
    pub value: f64,
}

impl std::fmt::Display for Limit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Limit { axis, op, value } = self;
        write!(f, "{axis} {op} {value}")
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Quantity {
    /// Copy number.
    Number(u64),
    /// Molarity (mol/L).
    Concentration(f64),
}

impl std::fmt::Display for Quantity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Number(count) => count.fmt(f),
            Self::Concentration(conc) => write!(f, "{conc}M"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub name: String,
    pub tag: Option<String>, // Should be 5-character ArrayString.
    pub quantity: Quantity,
    pub path: PathBuf,
    pub compartment_ids: Box<[String]>,
    pub rules: Option<Box<[String]>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Compartment {
    pub id: String,
    pub mask: Mask,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Mask {
    All,
    Voxels(PathBuf),
    Shape(Shape),
    Combination(Expr<String>),
}

#[derive(Debug, Clone, PartialEq)]
// TODO: This name is not appropriate.
pub enum Shape {
    Sphere { center: Center, radius: f32 }, // Consider the f64 situation.
    Cuboid { start: Anchor, end: Anchor },
}

impl Shape {
    fn space_filling_cuboid() -> Self {
        Self::Cuboid {
            start: Anchor::Start,
            end: Anchor::End,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Center {
    Center,
    Point(Point),
}

impl std::fmt::Display for Center {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Center => "center".fmt(f),
            Self::Point([x, y, z]) => write!(f, "{x}, {y}, {z}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Anchor {
    Start,
    Center,
    End,
    Point(Point),
}

impl std::fmt::Display for Anchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Start => "start".fmt(f),
            Self::Center => "center".fmt(f),
            Self::End => "end".fmt(f),
            Self::Point([x, y, z]) => write!(f, "{x}, {y}, {z}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Constraint {
    pub id: String,
    pub rule: Rule,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Rule {
    Limits(Expr<Limit>),
    Within { distance: f32, id: String },
    RotationAxes(Axes),
    Combination(Expr<String>),
}

#[derive(Debug, PartialEq)]
pub struct Config {
    pub general: General,
    pub space: Space,
    pub includes: Vec<PathBuf>,
    pub constraints: Vec<Constraint>,
    pub compartments: Vec<Compartment>,
    pub segments: Vec<Segment>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

impl Axes {
    pub fn list(&self) -> Box<[Axis]> {
        let mut v = Vec::with_capacity(3);
        if self.x {
            v.push(Axis::X);
        }
        if self.y {
            v.push(Axis::Y);
        }
        if self.z {
            v.push(Axis::Z);
        }
        v.into_boxed_slice()
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
