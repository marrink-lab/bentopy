use bentopy::core::config::legacy::CombinationExpression;

use crate::{mask::Mask, state::Compartment};

/// Apply some operation over the mask resulting from the first [`CombinationExpression`] in
/// `exprs` for each subsequent `expr`.
///
/// # Invariants
///
/// Assumes `exprs` has at least one entry.
fn fold<F: Fn(&mut Mask, Mask)>(
    exprs: &[CombinationExpression],
    compartments: &[Compartment],
    op: F,
) -> anyhow::Result<Mask> {
    // Invariant: Assumes `exprs` has at least one entry.
    let (first, exprs) = exprs.split_first().unwrap();
    let mut mask = execute(first, compartments)?;
    // Now apply `op` between the first mask and those resulting from the subsequent `exprs`.
    for expr in exprs {
        let src = execute(expr, compartments)?;
        op(&mut mask, src);
    }

    Ok(mask)
}

pub fn execute(expr: &CombinationExpression, compartments: &[Compartment]) -> anyhow::Result<Mask> {
    let output = match expr {
        // TODO: This is where a Cow may be great!
        CombinationExpression::Id(id) => {
            let compartment = compartments
                .iter()
                .find(|c| &c.id == id)
                .ok_or(anyhow::anyhow!("mask with id {id:?} not (yet) defined"))?;
            compartment.mask.clone()
        }
        CombinationExpression::Not(expr) => !execute(expr, compartments)?,
        CombinationExpression::Union(exprs) => fold(exprs, compartments, |mask, src| *mask &= src)?,
        CombinationExpression::Intersect(exprs) => {
            fold(exprs, compartments, |mask, src| *mask |= src)?
        }
    };

    Ok(output)
}
