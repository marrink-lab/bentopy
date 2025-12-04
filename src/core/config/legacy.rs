use std::{path::PathBuf, str::FromStr};

use serde::Deserialize;

pub use super::compartment_combinations::Expression as CombinationExpression;
use crate::core::config::{Axes, CompartmentID, Dimensions, defaults};

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
    defaults::BEAD_RADIUS as f32 // nm
}

fn max_tries_mult_default() -> u64 {
    defaults::MAX_TRIES_MULT
}

fn max_tries_rot_div_default() -> u64 {
    defaults::MAX_TRIES_ROT_DIV
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
    use std::collections::HashSet;

    use super::*;
    use crate::core::config::{self, legacy::convert::rule::parse_rule};

    impl Into<config::Expr<CompartmentID>> for CombinationExpression {
        fn into(self) -> config::Expr<CompartmentID> {
            type Expr = config::Expr<CompartmentID>;
            fn unflatten(
                expressions: Vec<CombinationExpression>,
                binary: impl Fn(Box<Expr>, Box<Expr>) -> Expr,
            ) -> Expr {
                // We recursively convert the children first.
                let exprs: Vec<config::Expr<CompartmentID>> = expressions
                    .into_iter()
                    .map(|expression| expression.into())
                    .collect();
                // Then, we stitch them together into a tree.
                let flat = exprs;
                flat.into_iter()
                    .reduce(|acc, item| binary(Box::new(acc), Box::new(item)))
                    .expect("TODO")
            }

            match self {
                Self::Id(id) => config::Expr::Term(id),
                Self::Not(expression) => config::Expr::Not(Box::new((*expression).into())),
                Self::Union(expressions) => unflatten(expressions, config::Expr::Or),
                Self::Intersect(expressions) => unflatten(expressions, config::Expr::And),
            }
        }
    }

    impl Into<config::Mask> for Mask {
        fn into(self) -> config::Mask {
            match self {
                Mask::Shape(shape) => config::Mask::Shape(match shape {
                    // TODO: This sucks but is true.
                    Shape::Spherical => panic!("a sphere without a radius is an undefined shape"),
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
                    radius: radius.expect("TODO"),
                }),
                Mask::Analytical {
                    shape: Shape::Cuboid | Shape::None,
                    center: _, // TODO: Error if not None.
                    radius: _, // TODO: Error if not None.
                } => config::Mask::Shape(config::Shape::space_filling_cuboid()),
                Mask::Voxels { path } => config::Mask::Voxels(path),
                Mask::Combination(expression) => {
                    // TODO: This conversion sucks and will be changed.
                    config::Mask::Combination(expression.into())
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
            // TODO: Implement parsing for rotation_axes properties.
            if rotation_axes == Default::default() {
                todo!("implement segment rotation axes for bent")
            }
            if initial_rotation != <[f32; 3]>::default() {
                unimplemented!("segment initial rotation is deprecated")
            }

            let rules = if rules.is_empty() {
                None
            } else {
                // This is a bit complicated. We need to
                // (1) come up with a label for each unique RuleExpression in the file,
                // (2) those rule expressions need to be formulated as a constraint and
                //     pushed to the constraints list. This is done in the conversion of
                //     legacy::Config to config::Config.
                let ids = rules
                    .iter()
                    .map(|re| {
                        let legacy_rule = parse_rule(re).expect("TODO");
                        rule::canonical_id(&legacy_rule)
                    })
                    // Only retain unique rules.
                    .collect::<HashSet<_>>();
                Some(ids.into_iter().collect())
            };

            config::Segment {
                name,
                tag,
                quantity: quantity.into(),
                path,
                compartment_ids: compartments.into_boxed_slice(),
                rules,
            }
        }
    }

    impl Into<config::Quantity> for Quantity {
        fn into(self) -> config::Quantity {
            match self {
                Self::Number(n) => config::Quantity::Number(n as u64),
                Self::Concentration(c) => config::Quantity::Concentration(c),
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
                        size,
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

            let constraints = segments
                .iter()
                .flat_map(|s| &s.rules)
                .map(|re| {
                    // First, we do a legacy parse of the rule.
                    let legacy_rule = parse_rule(re).expect("TODO");
                    // Come up with a unique id for this rule that will match how a segment converts
                    // the rule into the same id.
                    let id = rule::canonical_id(&legacy_rule);
                    // And convert the RuleExpression into a Rule.
                    let rule = legacy_rule.into();
                    config::Constraint { id, rule }
                })
                .collect();

            config::Config {
                general: config::General {
                    title: if title.is_empty() { None } else { Some(title) },
                    seed,
                    bead_radius: Some(bead_radius as f64),
                    max_tries_mult: Some(max_tries_mult),
                    max_tries_rot_div: Some(max_tries_rot_div),
                    rearrange_method: None,
                },
                space: config::Space {
                    dimensions: Some(size),
                    resolution: Some(resolution as f64),
                    periodic: Some(periodic),
                },
                includes: topol_includes
                    .unwrap_or_default()
                    .into_iter()
                    .map(Into::into)
                    .collect(),
                constraints,
                compartments: compartments.into_iter().map(Into::into).collect(),
                segments: segments.into_iter().map(Into::into).collect(),
            }
        }
    }

    // We're vendoring the Rule stuff until that can all be refactored out.
    mod rule {
        use std::num::ParseFloatError;
        use std::str::FromStr;

        use crate::core::config::{self, Axis, Limit, Op};

        use super::{CompartmentID, RuleExpression};

        // TODO: This should be part of the config parsing.
        pub fn parse_rule(expr: &RuleExpression) -> Result<Rule, ParseRuleError> {
            match expr {
                RuleExpression::Rule(s) => Rule::from_str(s),
                RuleExpression::Or(exprs) => Ok(Rule::Or(
                    exprs.iter().map(parse_rule).collect::<Result<_, _>>()?,
                )),
            }
        }

        pub fn canonical_id(rule: &Rule) -> String {
            match rule {
                Rule::Position(Limit { axis, op, value }) => {
                    let op = match op {
                        Op::SmallerThan => "lt",
                        Op::GreaterThan => "gt",
                    };
                    format!("{{lim/{axis:?}'{op}'{value}}}")
                }
                Rule::IsCloser(id, distance) => {
                    format!("{{win/{id}'{distance}}}")
                }
                Rule::Or(rules) => {
                    let ids = rules
                        .iter()
                        .map(|r| canonical_id(r))
                        .collect::<Vec<_>>()
                        .join("'");
                    format!("{{or/{ids}}}")
                }
            }
        }

        impl Into<config::Rule> for Rule {
            fn into(self) -> config::Rule {
                match self {
                    Rule::Position(limit) => config::Rule::Limits(config::Expr::Term(limit)),
                    Rule::IsCloser(id, distance) => config::Rule::Within { distance, id },
                    Rule::Or(_rules) => {
                        // TODO: This is quite easy to implement. Along the lines of
                        // compartment combinations, but just for rule ids.
                        todo!("rule combinations have not been implemented")
                    }
                }
            }
        }

        #[derive(Debug, Clone, PartialEq)]
        pub enum Rule {
            Position(Limit),
            IsCloser(CompartmentID, f32),

            /// A set of rules where any of them can be true for this [`Rule`] to apply.
            Or(Vec<Rule>),
        }

        impl FromStr for Rule {
            type Err = ParseRuleError;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let trimmed = s.trim();
                let mut words = trimmed.split_whitespace();
                let keyword = words.next().ok_or(ParseRuleError::Empty)?;
                match keyword {
                    kind @ ("less_than" | "greater_than") => {
                        let axis = words
                            .next()
                            .ok_or(ParseRuleError::SyntaxError("expected axis".to_string()))?
                            .parse()
                            .map_err(ParseRuleError::ParseAxisError)?;
                        let value = words
                            .next()
                            .ok_or(ParseRuleError::SyntaxError(
                                "expected scalar value".to_string(),
                            ))?
                            .parse()
                            .map_err(ParseRuleError::ParseScalarError)?;

                        let poscon = match kind {
                            "greater_than" => Limit {
                                axis,
                                op: Op::GreaterThan,
                                value,
                            },
                            "less_than" => Limit {
                                axis,
                                op: Op::SmallerThan,
                                value,
                            },
                            _ => unreachable!(), // By virtue of this branch's pattern.
                        };
                        Ok(Rule::Position(poscon))
                    }
                    "is_closer_to" => {
                        let compartment_id = words.next().ok_or(ParseRuleError::SyntaxError(
                            "expected compartment id".to_string(),
                        ))?;
                        let distance = words
                            .next()
                            .ok_or(ParseRuleError::SyntaxError(
                                "expected scalar value".to_string(),
                            ))?
                            .parse()
                            .map_err(ParseRuleError::ParseScalarError)?;

                        Ok(Rule::IsCloser(compartment_id.to_string(), distance))
                    }
                    unknown => Err(ParseRuleError::UnknownKeyword(unknown.to_string())),
                }
            }
        }

        #[derive(Debug, Clone)]
        pub enum ParseRuleError {
            Empty,
            UnknownKeyword(String),
            SyntaxError(String),
            ParseScalarError(ParseFloatError),
            ParseAxisError(String),
        }

        impl std::fmt::Display for ParseRuleError {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    ParseRuleError::Empty => write!(f, "no rule keyword was provided"),
                    ParseRuleError::UnknownKeyword(unknown) => {
                        write!(f, "encountered an unknown keyword: {unknown:?}")
                    }
                    ParseRuleError::SyntaxError(err) => write!(f, "syntax error: {err}"),
                    ParseRuleError::ParseScalarError(err) => {
                        write!(f, "could not parse float: {err}")
                    }
                    ParseRuleError::ParseAxisError(err) => write!(f, "could not parse axis: {err}"),
                }
            }
        }

        impl std::error::Error for ParseRuleError {}

        impl FromStr for Axis {
            type Err = String;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                match s {
                    "x" => Ok(Self::X),
                    "y" => Ok(Self::Y),
                    "z" => Ok(Self::Z),
                    weird => Err(format!(
                        "expected one of 'x', 'y', or 'z', but found {weird:?}"
                    )),
                }
            }
        }
    }
}
