//! Render a placement list to a gro file with speed.
//!
//! By Marieke Westendorp, 2024.
//! <ma3ke.cyber@gmail.com>
mod args;
mod limits;
mod render;
mod structure;

fn main() -> anyhow::Result<()> {
    render::render(clap::Parser::parse())
}
