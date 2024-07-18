use std::{io, path::Path};

use eightyseven::structure::BoxVecs;
use eightyseven::writer::WriteGro;
use glam::Vec3;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;

const WRITE_CHUNK_SIZE: usize = 0x1000000;
const FORMAT_CHUNK_SIZE: usize = WRITE_CHUNK_SIZE / 0x10;

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

pub trait WriteParExt {
    fn save_gro_par<P: AsRef<Path>>(&self, path: P) -> io::Result<()>;
}

impl WriteParExt for Structure {
    fn save_gro_par<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = io::BufWriter::new(file);
        write_par(self, &mut writer)
    }
}

pub fn write_par(structure: &Structure, writer: &mut impl io::Write) -> io::Result<()> {
    // Write the title.
    writeln!(writer, "solvated by {}", env!("CARGO_PKG_NAME"))?;

    // And the number of atoms in the system.
    let natoms = structure.natoms();
    writeln!(writer, "{natoms}")?;

    // Write the atoms.
    let mut n = 0;
    let start = std::time::Instant::now();
    for chunk in structure.atoms.chunks(WRITE_CHUNK_SIZE) {
        // Format this chunk in parallel.
        let formatted = chunk
            .par_chunks(FORMAT_CHUNK_SIZE)
            .map(|atoms| {
                atoms
                    .into_iter()
                    .map(Structure::format_atom_line)
                    .collect::<String>()
            })
            .collect_vec_list();
        // And write it away.
        for f in formatted.iter().flatten() {
            writer.write_all(f.as_bytes())?;
        }

        // Report the progress.
        n += chunk.len();
        let progress = n as f32 / natoms as f32;
        let percentage = progress * 100.0;
        let delta = std::time::Instant::now() - start;
        let time_left = delta.as_secs_f32() * (progress.recip() - 1.0);
        eprint!("\r({percentage:>4.1}%) Wrote {n:>9}/{natoms} atoms. (ETA {time_left:>4.0} s) ");
    }
    eprintln!();

    // End with the box vectors.
    // FIXME: What to do with >3 value box vectors?
    writeln!(writer, "{}", structure.boxvecs())?;

    writer.flush()
}
