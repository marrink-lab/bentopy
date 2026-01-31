use bentopy::core::config::{CompartmentID, Expr};

use crate::{mask::Mask, state::Compartment};

pub fn execute(expr: &Expr<CompartmentID>, compartments: &[Compartment]) -> anyhow::Result<Mask> {
    let output = match expr {
        Expr::Term(id) => {
            let compartment = compartments
                .iter()
                .find(|c| &c.id == id)
                .ok_or(anyhow::anyhow!("Mask with id {id:?} not (yet) defined"))?;
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
