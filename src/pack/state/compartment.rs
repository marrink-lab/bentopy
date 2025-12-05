use bentopy::core::config::{CompartmentID, Expr, Limit};

use crate::mask::{Dimensions, Mask};

pub struct Compartment {
    pub id: CompartmentID,
    pub mask: Mask,
}

pub fn distill_limits(expr: &Expr<Limit>, dimensions: Dimensions, resolution: f64) -> Mask {
    fn apply_limit(mut mask: Mask, limit: &Limit, resolution: f64) -> Mask {
        let &Limit { axis, op, value } = limit;
        match op {
            bentopy::core::config::Op::SmallerThan => {
                mask.apply_function(|pos| value < pos[axis as usize] as f64 * resolution)
            }
            bentopy::core::config::Op::GreaterThan => {
                mask.apply_function(|pos| value > pos[axis as usize] as f64 * resolution)
            }
        }

        mask
    }

    // TODO: This is naive. Also wrong perhaps. Let's test this thuroughly.
    fn d(expr: &Expr<Limit>, dimensions: Dimensions, resolution: f64) -> Mask {
        match expr {
            Expr::Term(limit) => apply_limit(Mask::new(dimensions), limit, resolution),
            Expr::Not(expr) => !d(expr, dimensions, resolution),
            Expr::Or(lhs, rhs) => {
                let mut m = d(lhs, dimensions, resolution);
                m &= d(rhs, dimensions, resolution);
                m
            }
            Expr::And(lhs, rhs) => {
                let mut m = d(lhs, dimensions, resolution);
                m |= d(rhs, dimensions, resolution);
                m
            }
        }
    }

    d(expr, dimensions, resolution)
}
