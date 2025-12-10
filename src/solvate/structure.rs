use std::io;

use eightyseven::structure::BoxVecs;
use eightyseven::writer::WriteGro;
use glam::Vec3;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{placement::PlaceMap, substitute::Substitution, water::WaterType};

pub type Structure = eightyseven::structure::Structure;

pub trait BoxVecsExtension {
    fn as_vec3(&self) -> glam::Vec3;
}

impl BoxVecsExtension for BoxVecs {
    fn as_vec3(&self) -> Vec3 {
        match self {
            // FIX: When a gro file with cuboid boundaries represented as a nine-valued vector
            // (trailing six values are zero) is loaded, we can canonicalize it as a three-valued
            // cuboid. This should probably be fixed in eightyseven, but I'll want to reconsider
            // this and related issues properly at another moment.
            &BoxVecs::Short(three) | &BoxVecs::Full([three @ .., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) => {
                Vec3::from_array(three)
            }
            BoxVecs::Full(nine) => todo!("nine-value boxvecs ({nine:?}) are currently unsupported"),
        }
    }
}

pub fn write_structure<const PAR: bool>(
    writer: &mut impl io::Write,
    structure: &Structure,
    solvent_placemap: &PlaceMap,
    substitutions: &[Substitution],
    water_type: WaterType,
    buffer_size: usize,
) -> io::Result<()> {
    let natoms_structure = structure.natoms();
    let natoms_solvent = solvent_placemap.unoccupied_count() as usize * water_type.total_natoms();
    let natoms_substitutes = substitutions.iter().map(|s| s.natoms()).sum::<usize>();
    let natoms = natoms_structure + natoms_solvent + natoms_substitutes;

    // Write the gro header.
    writeln!(
        writer,
        "solvated by {} (v{})",
        env!("CARGO_BIN_NAME"),
        bentopy::core::version::VERSION
    )?;
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
    let mut placements = solvent_placemap.iter_atoms_chunks();
    let mut buffer = Vec::new();
    let mut n = 0;
    'write_atoms: loop {
        buffer.clear();

        'fill_buffer: while buffer.len() < buffer_size {
            match placements.next() {
                Some(sp) => buffer.extend(sp),
                None if buffer.is_empty() => break 'write_atoms,
                None => break 'fill_buffer,
            }
        }

        // TODO: Consider breaking this into specific implementations for vecs of waters and single
        // waters. Don't want the performance of the usual single-bead martini water to suffer!
        let beads = &buffer;
        if PAR {
            let lines: String = beads
                .par_iter()
                .flat_map(|bead| water_type.expand(bead.position))
                .map(|a| Structure::format_atom_line(&a))
                .collect();
            writer.write_all(lines.as_bytes())?;
        } else {
            let expanded = beads
                .iter()
                .flat_map(|bead| water_type.expand(bead.position));
            for atom in expanded {
                let line = Structure::format_atom_line(&atom);
                writer.write_all(line.as_bytes())?;
            }
        }

        // Report the progress.
        n += buffer.len() * water_type.total_natoms();
        let progress = n as f32 / natoms_solvent as f32;
        let percentage = progress * 100.0;
        let delta = std::time::Instant::now() - start;
        let time_left = delta.as_secs_f32() * (progress.recip() - 1.0);
        eprint!(
            "\r({percentage:>4.1}%) Wrote {n:>9}/{natoms_solvent} solvent atoms. (ETA {time_left:>4.0} s) "
        );
    }
    eprintln!();
    eprintln!(
        "Done writing solvent atoms, took {:.3} s.",
        start.elapsed().as_secs_f32()
    );

    if !substitutions.is_empty() {
        // Write out the substitutes at the end of the structure.
        let start = std::time::Instant::now();
        eprintln!("Writing substitutes        ({natoms_substitutes:>9}/{natoms} total atoms)...");
        for substitution in substitutions {
            let name = substitution.name();
            let number = substitution.natoms();
            eprint!("\t{name}, {number} atoms... ");
            let start = std::time::Instant::now();
            let mut placements = substitution.iter_atoms();
            loop {
                buffer.clear();
                buffer.extend(placements.by_ref().take(buffer_size));
                if buffer.is_empty() {
                    break;
                }

                let beads = &buffer;
                if PAR {
                    let lines: String = beads.par_iter().map(Structure::format_atom_line).collect();
                    writer.write_all(lines.as_bytes())?;
                } else {
                    for atom in beads {
                        let line = Structure::format_atom_line(atom);
                        writer.write_all(line.as_bytes())?;
                    }
                }
            }
            eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
        }
        eprintln!(
            "Done writing substitutes, took {:.3} s.",
            start.elapsed().as_secs_f32()
        );
    }

    // Finally, we write the box vectors.
    writeln!(writer, "{}", structure.boxvecs())?;

    writer.flush()?;

    Ok(())
}
