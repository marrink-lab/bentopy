use anyhow::Result;

use std::path::PathBuf;

pub mod bent;
mod compartment_combinations;
pub mod legacy;

pub mod defaults {
    pub const TITLE: &str = "Bentopy system";
    /// Bead radius for voxelization in nm.
    pub const BEAD_RADIUS: f64 = 0.20;
    pub const PERIODIC: bool = true;
    pub const MAX_TRIES_MULT: u64 = 1000;
    pub const MAX_TRIES_ROT_DIV: u64 = 100;
}

/// Avogadro's number (per mol).
const N_A: f64 = 6.0221415e23;

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
            _ => Err("unknown axis"),
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
    LessThan,
    GreaterThan,
}

impl Op {
    fn reverse(self) -> Self {
        match self {
            Self::LessThan => Self::GreaterThan,
            Self::GreaterThan => Self::LessThan,
        }
    }
}

impl TryFrom<char> for Op {
    type Error = &'static str;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '<' => Ok(Self::LessThan),
            '>' => Ok(Self::GreaterThan),
            _ => Err("unknown operator"),
        }
    }
}

impl std::fmt::Display for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LessThan => '<',
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

impl Quantity {
    /// Determine the number of segments that is implied by this [`Quantity`].
    ///
    /// In case this `Quantity` is a [`Quantity::Concentration`], the number of segments is
    /// lazily determined from the provided `volume`, and rounded.
    ///
    /// The value returned by `volume` must be in cubic nanometers (nm³).
    pub fn bake<F: Fn() -> f64>(&self, volume: F) -> u64 {
        match *self {
            Quantity::Number(n) => n,
            Quantity::Concentration(c) => {
                let v = volume() * 1e-24; // From nm³ to L.
                let n = N_A * c * v;
                f64::round(n) as u64
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
    pub rules: Box<[String]>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Compartment {
    pub id: CompartmentID,
    pub mask: Mask,
}

impl Compartment {
    /// Returns whether this [`Compartment`] is dependent on other compartments.
    pub fn is_predefined(&self) -> bool {
        match &self.mask {
            Mask::All | Mask::Voxels(_) | Mask::Shape(_) | Mask::Limits(_) => true,
            Mask::Within { .. } | Mask::Combination(_) => false,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Mask {
    All,
    Voxels(PathBuf),
    Shape(Shape),
    Limits(Expr<Limit>),
    // These are constructed by referencing previously defined masks.
    Within { distance: f32, id: CompartmentID },
    Combination(Expr<CompartmentID>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Shape {
    Sphere { center: Center, radius: f32 }, // Consider the f64 situation.
    Cuboid { start: Anchor, end: Anchor },
    // TODO: More?
}

impl std::fmt::Display for Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Sphere { center, radius } => {
                write!(f, "sphere at {center} with radius {radius}")
            }
            Self::Cuboid { start, end } => write!(f, "cuboid from {start} to {end}"),
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Constraint {
    pub id: String,
    pub rule: Rule,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Rule {
    RotationAxes(Axes),
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
