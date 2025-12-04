use std::num::ParseFloatError;
use std::str::FromStr;

use bentopy::core::config::legacy::RuleExpression;

use crate::mask::{Dimensions, Mask};
use crate::state::{Compartment, CompartmentID};

type Scalar = f32;

// TODO: This should be part of the config parsing.
pub fn parse_rule(expr: &RuleExpression) -> Result<Rule, ParseRuleError> {
    match expr {
        RuleExpression::Rule(s) => Rule::from_str(s),
        RuleExpression::Or(exprs) => Ok(Rule::Or(
            exprs.iter().map(parse_rule).collect::<Result<_, _>>()?,
        )),
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Rule {
    Position(PositionConstraint),
    IsCloser(CloseStyleBikeshed, CompartmentID, Scalar),

    /// A set of rules where any of them can be true for this [`Rule`] to apply.
    Or(Vec<Rule>),
}

impl Rule {
    fn distill(
        &self,
        dimensions: Dimensions,
        resolution: f32,
        compartments: &[Compartment],
    ) -> Mask {
        match self {
            Rule::Position(poscon) => {
                let mut mask = Mask::new(dimensions);
                let axi = poscon.axis().as_idx();
                match poscon {
                    &PositionConstraint::GreaterThan(_, value) => {
                        mask.apply_function(|pos| value < pos[axi] as f32 * resolution)
                    }
                    &PositionConstraint::LessThan(_, value) => {
                        mask.apply_function(|pos| value > pos[axi] as f32 * resolution)
                    }
                }

                mask
            }
            Rule::IsCloser(_style, compartment_id, distance) => {
                eprintln!("NOTE: Doing the whole mask thing right now.");
                let compartment = compartments
                    .iter()
                    .find(|c| c.id == compartment_id.as_str())
                    // FIXME: This condition can be checked at time of State creation.
                    .expect("the compartment ID in the rule must exist");

                let m = compartment.get_distance_mask(distance / resolution);
                eprintln!("NOTE: Done.");
                m
            }
            Rule::Or(rules) => {
                let mut mask = Mask::new(dimensions);
                // FIXME: I'd rather just apply it to the same mask once.
                for rule in rules {
                    mask |= rule.distill(dimensions, resolution, compartments)
                }

                mask
            }
        }
    }
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
                    "less_than" => PositionConstraint::LessThan(axis, value),
                    "greater_than" => PositionConstraint::GreaterThan(axis, value),
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

                Ok(Rule::IsCloser(
                    CloseStyleBikeshed::BoxCenter,
                    compartment_id.to_string(),
                    distance,
                ))
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
            ParseRuleError::ParseScalarError(err) => write!(f, "could not parse float: {err}"),
            ParseRuleError::ParseAxisError(err) => write!(f, "could not parse axis: {err}"),
        }
    }
}

impl std::error::Error for ParseRuleError {}

pub fn distill(
    rules: &[Rule],
    dimensions: Dimensions,
    resolution: f32,
    compartments: &[Compartment],
) -> Mask {
    let mut mask = Mask::fill::<true>(dimensions);
    for rule in rules {
        mask &= rule.distill(dimensions, resolution, compartments);
    }

    mask
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X = 0,
    Y = 1,
    Z = 2,
}

impl Axis {
    pub const fn as_idx(&self) -> usize {
        *self as usize
    }
}

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

#[derive(Debug, Clone, PartialEq)]
pub enum PositionConstraint {
    GreaterThan(Axis, Scalar),
    LessThan(Axis, Scalar),
}

impl PositionConstraint {
    pub fn axis(&self) -> Axis {
        match self {
            PositionConstraint::GreaterThan(axis, _) | PositionConstraint::LessThan(axis, _) => {
                *axis
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CloseStyleBikeshed {
    BoxCenter,
    // AnyFilledVoxel,
}
