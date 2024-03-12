#![allow(dead_code)]

use std::{path::PathBuf, str::FromStr};

use pyo3::prelude::*;
use render_placements::{render, Limits, Mode};

#[pyfunction]
fn py_render_placements(
    input_path: PathBuf,
    output_path: PathBuf,
    topol_path: Option<PathBuf>,
    root: Option<PathBuf>,
    limits: Option<String>,
    mode: Option<String>,
) -> std::io::Result<()> {
    eprintln!("Hello, world");
    let mode = mode
        .and_then(|m| Mode::from_str(&m).ok())
        .unwrap_or_default();
    let limits = limits.and_then(|l| Limits::from_str(&l).ok());
    render(input_path, output_path, topol_path, root, limits, mode)
}

/// Read xtc files, fast.
#[pymodule]
fn _render(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_render_placements, m)?)?;

    Ok(())
}
