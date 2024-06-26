use glam::Vec3;

#[derive(Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z,
}

type Scalar = f32;

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

pub enum Rule {
    Position(PositionConstraint),
}

impl Rule {
    /// Return whether this rule is satisfied or not.
    ///
    /// The provided `pos` must be in nanometers.
    pub fn is_satisfied(&self, min: Vec3, max: Vec3) -> bool {
        match self {
            Rule::Position(poscon) => {
                // TODO: Gotta find a way to know the bounding box or all positions
                // of the structure placement. We calculate these somewhere already
                // maybe??
                let (min, max) = match poscon.axis() {
                    Axis::X => (min.x, max.x),
                    Axis::Y => (min.y, max.y),
                    Axis::Z => (min.z, max.z),
                };
                match poscon {
                    &PositionConstraint::GreaterThan(_, value) => value > max,
                    &PositionConstraint::LessThan(_, value) => value < min,
                }
            }
        }
    }

    fn is_lightweight(&self) -> bool {
        match self {
            Rule::Position(_) => true,
        }
    }
}

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
