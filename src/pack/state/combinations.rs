use bentopy::core::config::{CompartmentID, Expr, legacy::CombinationExpression};

use crate::{mask::Mask, state::Compartment};

/// Apply some operation over the mask resulting from the first [`CombinationExpression`] in
/// `exprs` for each subsequent `expr`.
///
/// # Invariants
///
/// Assumes `exprs` has at least one entry.
fn legacy_fold<F: Fn(&mut Mask, Mask)>(
    exprs: &[CombinationExpression],
    compartments: &[Compartment],
    op: F,
) -> anyhow::Result<Mask> {
    // Invariant: Assumes `exprs` has at least one entry.
    let (first, exprs) = exprs.split_first().unwrap();
    let mut mask = legacy_execute(first, compartments)?;
    // Now apply `op` between the first mask and those resulting from the subsequent `exprs`.
    for expr in exprs {
        let src = legacy_execute(expr, compartments)?;
        op(&mut mask, src);
    }

    Ok(mask)
}

pub fn legacy_execute(
    expr: &CombinationExpression,
    compartments: &[Compartment],
) -> anyhow::Result<Mask> {
    let output = match expr {
        // TODO: This is where a Cow may be great!
        CombinationExpression::Id(id) => {
            let compartment = compartments
                .iter()
                .find(|c| &c.id == id)
                .ok_or(anyhow::anyhow!("mask with id {id:?} not (yet) defined"))?;
            compartment.mask.clone()
        }
        CombinationExpression::Not(expr) => !legacy_execute(expr, compartments)?,
        CombinationExpression::Union(exprs) => {
            legacy_fold(exprs, compartments, |mask, src| *mask &= src)?
        }
        CombinationExpression::Intersect(exprs) => {
            legacy_fold(exprs, compartments, |mask, src| *mask |= src)?
        }
    };

    Ok(output)
}

pub fn execute(expr: &Expr<CompartmentID>, compartments: &[Compartment]) -> anyhow::Result<Mask> {
    let output = match expr {
        Expr::Term(id) => {
            let compartment = compartments
                .iter()
                .find(|c| &c.id == id)
                .ok_or(anyhow::anyhow!("mask with id {id:?} not (yet) defined"))?;
            compartment.mask.clone()
        }
        Expr::Not(expr) => !execute(expr, compartments)?,
        // TODO: This is quite a horrible, no-good transliteration of the legacy_execute function.
        // There is some much here that can be improved upon, in terms of memory efficiency, for
        // example.
        Expr::Or(lhs, rhs) => {
            let mut m = execute(lhs, compartments)?;
            m &= execute(rhs, compartments)?;
            m
        }
        Expr::And(lhs, rhs) => {
            let mut m = execute(lhs, compartments)?;
            m |= execute(rhs, compartments)?;
            m
        }
    };

    Ok(output)
}
