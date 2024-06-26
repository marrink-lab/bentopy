use glam::Vec3;

#[derive(Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z,
}

pub enum PositionConstraint {
    GreaterThan(Axis, f64),
    LessThan(Axis, f64),
    At(Axis, f64),
}

impl PositionConstraint {
    pub fn axis(&self) -> Axis {
        match self {
            PositionConstraint::GreaterThan(axis, _)
            | PositionConstraint::LessThan(axis, _)
            | PositionConstraint::At(axis, _) => *axis,
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
    pub fn is_satisfied(&self, pos: Vec3) -> bool {
        match self {
            Rule::Position(poscon) => {
                // TODO: Gotta find a way to know the bounding box or all positions
                // of the structure placement. We calculate these somewhere already
                // maybe??
                let p = match poscon.axis() {
                    Axis::X => pos.x,
                    Axis::Y => pos.y,
                    Axis::Z => pos.z,
                } as f64;
                match poscon {
                    &PositionConstraint::GreaterThan(_, value) => value > p,
                    &PositionConstraint::LessThan(_, value) => value < p,
                    // FIXME: Have a (resolution-sized?) grace edge for the At.
                    &PositionConstraint::At(_, value) => value == p,
                }
            }
        }
    }
}
