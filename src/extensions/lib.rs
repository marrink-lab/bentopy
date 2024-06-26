use pyo3::prelude::*;

mod place;
mod render;
mod voxelize;

#[pymodule]
fn _extensions(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(place::py_place, m)?)?;
    m.add_function(wrap_pyfunction!(render::py_render_placements, m)?)?;
    m.add_function(wrap_pyfunction!(voxelize::py_voxelize, m)?)?;

    Ok(())
}
