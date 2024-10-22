use std::num::ParseFloatError;
use std::ops::{BitAndAssign, BitOrAssign, Not};
use std::str::FromStr;

use glam::U64Vec3;

use crate::mask::{Dimensions, Mask, Position};
use crate::state::{Compartment, CompartmentID};

type Scalar = f32;

#[derive(Debug, Clone, PartialEq)]
pub enum Rule {
    Position(PositionConstraint),
    IsCloser(CloseStyleBikeshed, CompartmentID, Scalar),

    /// A set of rules where any of them can be true for this [`Rule`] to apply.
    Or(Vec<Rule>),
}

impl Rule {
    /// Return whether this rule is satisfied or not.
    ///
    /// The provided `pos` must be in nanometers.
    pub fn is_satisfied(
        &self,
        position: Position,
        resolution: f32,
        voxels: &Mask,
        compartments: &[Compartment],
    ) -> bool {
        match self {
            Rule::Position(poscon) => {
                let min = U64Vec3::from_array(position).as_vec3() * resolution;
                let max = min + U64Vec3::from_array(voxels.dimensions()).as_vec3() * resolution;
                let (min, max) = match poscon.axis() {
                    Axis::X => (min.x, max.x),
                    Axis::Y => (min.y, max.y),
                    Axis::Z => (min.z, max.z),
                };
                match poscon {
                    &PositionConstraint::GreaterThan(_, value) => value < min,
                    &PositionConstraint::LessThan(_, value) => value > max,
                }
            }
            Rule::IsCloser(CloseStyleBikeshed::BoxCenter, compartment_id, distance) => {
                assert!(*distance > 0.0, "distance must be a positive float");

                // FIXME: Compartments really ought to be a HashMap of IDs to Compartments.
                let compartment = compartments
                    .iter()
                    .find(|c| c.id == compartment_id.as_str())
                    // FIXME: This condition can be checked at time of State creation.
                    .expect("the compartment ID in the rule must exist");

                let voxel_radius = distance / resolution;
                let distance_mask = compartment.get_distance_mask(voxel_radius);

                // FIXME: Basically the internals of Session::check_collisions here.
                let occupied_indices = voxels.indices_where::<true>().map(U64Vec3::from_array);
                let position = U64Vec3::from_array(position);
                let voxels_dimensions = U64Vec3::from_array(voxels.dimensions());
                let voxels_center = position + voxels_dimensions / 2;
                occupied_indices
                    .map(|p| (p + voxels_center).to_array())
                    .all(|idx| !distance_mask.get_periodic(idx))
            }
            Rule::Or(rules) => rules
                .iter()
                .any(|rule| rule.is_satisfied(position, resolution, voxels, compartments)),
        }
    }

    fn is_lightweight(&self) -> bool {
        match self {
            Rule::Position(_) => true,
            Rule::IsCloser(_, _, _) => false,
            Rule::Or(rules) => rules.iter().all(Rule::is_lightweight),
        }
    }

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

impl BitOrAssign for Mask {
    fn bitor_assign(&mut self, rhs: Self) {
        assert_eq!(
            self.dimensions(),
            rhs.dimensions(),
            "the dimensions of both masks must be identical"
        );
        // For good measure, so the compiler gets it.
        assert_eq!(self.n_backings(), rhs.n_backings()); // FIXME: Is this one necessary?

        self.backings
            .iter_mut()
            .zip(rhs.backings.iter())
            .for_each(|(s, &m)| *s |= m);
    }
}

impl BitAndAssign for Mask {
    // TODO: Perhaps introduce a macro to set up functions like this?
    fn bitand_assign(&mut self, rhs: Self) {
        assert_eq!(
            self.dimensions(),
            rhs.dimensions(),
            "the dimensions of both masks must be identical"
        );
        // For good measure, so the compiler gets it.
        assert_eq!(self.n_backings(), rhs.n_backings()); // FIXME: Is this one necessary?

        self.backings
            .iter_mut()
            .zip(rhs.backings.iter())
            .for_each(|(s, &m)| *s &= m);
    }
}

impl Not for Mask {
    type Output = Self;

    fn not(mut self) -> Self::Output {
        self.backings.iter_mut().for_each(|b| *b = !*b);
        self
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
        write!(f, "{self:?}") // FIXME: This can made much more pretty and clear.
    }
}

impl std::error::Error for ParseRuleError {}

/// Split rules into a lightweight and heavier subset.
///
/// Split a set of rules into those that are cheap to compute for some position, and those that are
/// more expensive and that should be considered only for promising candidates, such as after
/// collision checking.
pub fn split<'r>(
    rules: &'r [Rule],
) -> (
    impl Iterator<Item = &'r Rule>,
    impl Iterator<Item = &'r Rule>,
) {
    fn criterion(rule: &Rule) -> bool {
        rule.is_lightweight()
    }
    (
        rules.iter().filter(|rule| criterion(rule)),
        rules.iter().filter(|rule| !criterion(rule)),
    )
}

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
