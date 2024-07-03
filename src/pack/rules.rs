use std::num::ParseFloatError;
use std::str::FromStr;

use glam::U64Vec3;

use crate::mask::{Mask, Position};
use crate::state::{Compartment, CompartmentID};

type Scalar = f32;

#[derive(Clone)]
pub enum Rule {
    Position(PositionConstraint),
    IsCloser(CloseStyleBikeshed, CompartmentID, Scalar),
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

                let position = U64Vec3::from_array(position);
                let voxels_dimensions = U64Vec3::from_array(voxels.dimensions());
                let voxels_center = position + voxels_dimensions / 2;

                let voxel_radius = distance / resolution;

                // TODO: All of the asserts littered here must be dealt with earlier on. Sometimes
                // they do describe artificial limitations, that just need some more implementation
                // of literal edge cases.

                assert!(voxel_radius.powi(2) as u64 <= voxels_dimensions.length_squared());
                let inscribed_box_size = U64Vec3::splat((voxel_radius * 2.0).ceil() as u64);
                let circumscribed_box_size =
                    U64Vec3::splat((voxel_radius * f32::sqrt(1.0 / 3.0)).ceil() as u64);
                let inscribed_box_center = inscribed_box_size / 2;
                let circumscribed_box_center = circumscribed_box_size / 2;
                // Check whether the upcoming subtractions are safe.
                assert!(inscribed_box_center.cmple(voxels_center).all());
                assert!(circumscribed_box_center.cmple(voxels_center).all());
                let inscribed_box_start = voxels_center - inscribed_box_center;
                let circumscribed_box_start = voxels_center - circumscribed_box_center;

                // FIXME: We're assuming that the box sizes don't go beyond the actual mask size
                // from some position. That's really bad.

                let inscribed = compartment.mask.slice_from_start_dimensions(
                    inscribed_box_start.to_array(),
                    inscribed_box_size.to_array(),
                );

                // Since this is the cube inscribed by the radius, any true voxel we find implies
                // that the distance criterion is met! So we are satisfied with just one hit :)
                if inscribed.any::<true>() {
                    // Found a hit! There is at least one voxel from the compartment mask that is
                    // at least `distance` from the `position`.
                    return true;
                }

                // TODO: We should skip checking over the sections we already checked.

                let position_f = position.as_vec3(); // FIXME: Horrible names.
                let voxel_radius2 = voxel_radius.powi(2);
                // Haven't found a hit yet, but maybe if we expand our search field a bit.
                let circumscribed = compartment.mask.slice_from_start_dimensions(
                    circumscribed_box_start.to_array(),
                    circumscribed_box_size.to_array(),
                );
                // In this case, we actually do need to check the distance, because in the corners
                // of this cube, the distance to our point of interest is actually greater than our
                // query `distance`.
                if circumscribed.iter_where::<true>().any(|idx| {
                    let hit = U64Vec3::from_array(idx).as_vec3();
                    let d2 = position_f.distance_squared(hit);
                    d2 <= voxel_radius2
                }) {
                    // Found a hit! There is at least one voxel from the compartment mask that is
                    // at least `distance` from the `position`.
                    return true;
                }

                false
            }
        }
    }

    fn is_lightweight(&self) -> bool {
        match self {
            Rule::Position(_) => true,
            Rule::IsCloser(_, _, _) => false,
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
                    _ => unreachable!(), // By virtue of this branche's pattern.
                };
                Ok(Rule::Position(poscon))
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
    impl Iterator<Item = &Rule> + 'r,
    impl Iterator<Item = &Rule> + 'r,
) {
    fn criterion(rule: &Rule) -> bool {
        rule.is_lightweight()
    }
    (
        rules.iter().filter(|rule| criterion(rule)),
        rules.iter().filter(|rule| !criterion(rule)),
    )
}

#[derive(Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z,
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

#[derive(Clone)]
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

#[derive(Clone)]
pub enum CloseStyleBikeshed {
    BoxCenter,
    // AnyFilledVoxel,
}
