//! Render a placement list to a gro file with speed.
//!
//! By Marieke Westendorp, 2024.
//! <ma3ke.cyber@gmail.com>
mod args;
mod limits;
mod render;
mod structure;

fn main() -> anyhow::Result<()> {
    let args::Args {
        input,
        output,
        topol,
        root,
        mode,
        limits,
        resnum_mode,
        ignore_tags,
        verbose,
    } = clap::Parser::parse();

    render::render(
        input,
        output,
        topol,
        root,
        limits,
        mode,
        resnum_mode,
        ignore_tags,
        verbose,
    )
}
