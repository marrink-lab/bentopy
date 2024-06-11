use std::{path::PathBuf, str::FromStr};

use pyo3::prelude::*;
use render_placements::{render, Limits, Mode, ResnumMode};

#[pyfunction]
pub(crate) fn py_render_placements(
    input_path: PathBuf,
    output_path: PathBuf,
    topol_path: Option<PathBuf>,
    root: Option<PathBuf>,
    limits: Option<String>,
    mode: Option<String>,
    resnum_mode: Option<String>,
) -> std::io::Result<()> {
    eprintln!("Hello, world");
    let mode = mode
        .and_then(|m| Mode::from_str(&m).ok())
        .unwrap_or_default();
    let resnum_mode = resnum_mode
        .and_then(|m| ResnumMode::from_str(&m).ok())
        .unwrap_or_default();
    let limits = limits.and_then(|l| Limits::from_str(&l).ok());
    render(
        input_path,
        output_path,
        topol_path,
        root,
        limits,
        mode,
        resnum_mode,
    )
}
