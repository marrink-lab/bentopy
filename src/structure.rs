use std::io;

use eightyseven::structure::BoxVecs;
use eightyseven::writer::WriteGro;
use glam::UVec3;
use glam::Vec3;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::placement::PlaceMap;

pub type Structure = eightyseven::structure::Structure;

pub trait BoxVecsExtension {
    fn as_vec3(&self) -> glam::Vec3;
}

impl BoxVecsExtension for BoxVecs {
    fn as_vec3(&self) -> Vec3 {
        match self {
            &BoxVecs::Short(three) => Vec3::from_array(three),
            BoxVecs::Full(nine) => todo!("nine-value boxvecs ({nine:?}) are currently unsupported"),
        }
    }
}

/// Iterate over all positions between `(0, 0, 0)` and `dimensions`.
fn iter_3d(dimensions: UVec3) -> impl Iterator<Item = UVec3> {
    let [dx, dy, dz] = dimensions.to_array();
    (0..dz).flat_map(move |z| (0..dy).flat_map(move |y| (0..dx).map(move |x| UVec3::new(x, y, z))))
}

pub fn write_structure<const PAR: bool>(
    writer: &mut impl io::Write,
    structure: &Structure,
    template: &Structure,
    placemap: &PlaceMap,
    buffer_size: usize,
) -> io::Result<()> {
    let natoms_structure = structure.natoms();
    let natoms_solvent = placemap.unoccupied_count();
    let natoms = natoms_structure + natoms_solvent;

    // Write the gro header.
    writeln!(writer, "solvated by {}", env!("CARGO_PKG_NAME"))?;
    writeln!(writer, "{natoms}")?;

    // First, we write the original structure.
    let start = std::time::Instant::now();
    eprint!("Writing original structure ({natoms_structure:>9}/{natoms} total atoms)... ");
    for line in structure.format_atom_lines_iter() {
        writer.write_all(line.as_bytes())?;
    }
    eprintln!("Done, took {:.3} s.", start.elapsed().as_secs_f32());

    // Now for the solvent. We will go over the cells described by the placemap, and format and
    // write the atom lines per cell. This means that we can avoid creating the whole solvent
    // structure at once, which would be tremendously memory intensive.
    let start = std::time::Instant::now();
    eprintln!("Writing solvent structure  ({natoms_solvent:>9}/{natoms} total atoms)...");
    let offset = |cell_pos: UVec3| cell_pos.as_vec3() * template.boxvecs.as_vec3();
    let dimensions = placemap.dimensions();
    let mut sps = iter_3d(dimensions).map(|cell_pos| {
        let translation = offset(cell_pos);
        let placement = placemap.get(cell_pos).unwrap();
        let solvent_positions =
            placemap
                .solvent
                .atoms()
                .zip(placement)
                .filter_map(move |(&sb, occupied)| {
                    if occupied {
                        None
                    } else {
                        let mut sb = sb; // Copy the solvent bead.
                        sb.position += translation;
                        Some(sb)
                    }
                });
        solvent_positions.collect::<Box<[_]>>()
    });

    let mut buffer = Vec::new();
    let mut n = 0;
    'write_atoms: loop {
        buffer.clear();

        'fill_buffer: while buffer.len() < buffer_size {
            match sps.next() {
                Some(sp) => buffer.extend(sp),
                None if buffer.len() == 0 => break 'write_atoms,
                None => break 'fill_buffer,
            }
        }

        let atoms = &buffer;
        if PAR {
            let lines: String = atoms.par_iter().map(Structure::format_atom_line).collect();
            writer.write_all(lines.as_bytes())?;
        } else {
            for atom in atoms {
                let line = Structure::format_atom_line(&atom);
                writer.write_all(line.as_bytes())?;
            }
        }

        // Report the progress.
        n += buffer.len();
        let progress = n as f32 / natoms_solvent as f32;
        let percentage = progress * 100.0;
        let delta = std::time::Instant::now() - start;
        let time_left = delta.as_secs_f32() * (progress.recip() - 1.0);
        eprint!("\r({percentage:>4.1}%) Wrote {n:>9}/{natoms_solvent} solvent atoms. (ETA {time_left:>4.0} s) ");
    }
    eprintln!();
    eprintln!(
        "Done writing solvent atoms, took {:.3} s.",
        start.elapsed().as_secs_f32()
    );

    // Eventually, this is where the ions would be written to the file. But we're not there yet.

    // Finally, we write the box vectors.
    writeln!(writer, "{}", structure.boxvecs())?;

    writer.flush()?;

    Ok(())
}
