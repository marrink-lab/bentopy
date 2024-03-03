use std::{
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
};

use arraystring::{typenum::U5, ArrayString};
use glam::Vec3;
use rayon::prelude::*;

/// The maximum number that can be represented by a 5-digit integer as found in gro files is 99999.
/// Any value up to this number can be correctly displayed within 5 characters.
const GRO_INTEGER_LIMIT: u32 = 100000;

type GroStr = ArrayString<U5>;

#[derive(Debug, Clone, Copy, Default)]
pub struct Bead {
    pub resnum: u32,
    pub resname: GroStr,
    pub atomname: GroStr,
    pub atomnum: u32,
    pub pos: Vec3,
}

impl Bead {
    pub fn new(
        resnum: u32,
        resname: &str,
        atomname: &str,
        atomnum: u32,
        pos: impl Into<Vec3>,
    ) -> Self {
        Self {
            resnum,
            resname: resname.into(),
            atomname: atomname.into(),
            atomnum,
            pos: pos.into(),
        }
    }

    // TODO: Consider whether inlining like this is actually appropriate?
    #[inline(always)]
    pub fn format_gro(&self, s: &mut String) {
        let resnum = self.resnum % GRO_INTEGER_LIMIT;
        let resname = self.resname.to_string(); // TODO: The ArrayString size problem remains.
        let atomname = self.atomname.to_string();
        let atomnum = self.atomnum % GRO_INTEGER_LIMIT;
        let [x, y, z] = self.pos.to_array();

        s.push_str(&format!(
            "{resnum:>5}{resname:<5}{atomname:>5}{atomnum:>5}{x:>8.3}{y:>8.3}{z:>8.3}"
        ));
    }

    // TODO: Consider whether inlining like this is actually appropriate?
    #[inline(always)]
    pub fn write_gro(&self, mut writer: impl Write) -> io::Result<()> {
        let resnum = self.resnum % GRO_INTEGER_LIMIT;
        let resname = self.resname.to_string(); // TODO: The ArrayString size problem remains.
        let atomname = self.atomname.to_string();
        let atomnum = self.atomnum % GRO_INTEGER_LIMIT;
        let [x, y, z] = self.pos.to_array();
        writeln!(
            writer,
            "{resnum:>5}{resname:<5}{atomname:>5}{atomnum:>5}{x:>8.3}{y:>8.3}{z:>8.3}"
        )
    }
}

#[derive(Debug, Default)]
pub struct Structure {
    pub beads: Vec<Bead>,
    pub boxvec: Vec3,
}

impl Structure {
    pub fn read_from_gro(reader: impl Read) -> io::Result<Self> {
        let reader = BufReader::new(reader);
        let mut lines = reader.lines();

        // Skip the title.
        let _title = lines.next();

        // Read the number of upcoming atoms.
        let n_atoms: usize = lines
            .next()
            .expect("unexpected end of file: no number of atoms specified")?
            .trim()
            .parse()
            .expect("invalid integer: could not parse number of atoms");

        // Stream in the firehose of atoms.
        let mut beads = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let line = lines
                .next()
                .expect("unexpected end of file: still reading atoms")?;

            let resnum = line[0..5].trim().parse().expect("invalid resnum integer");
            let resname = line[5..10].trim();
            let atomname = line[10..15].trim();
            let atomnum = line[15..20]
                .trim()
                .parse()
                .expect("invalid atomnum integer");
            // NOTE: This position parsing may fail when the last float does not fill out all its
            // decimal zeros. The gro file documentation explicitly describes the C format string
            // that describes an atom line ("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"). The C
            // behavior here is to print trailing zeros. I take that to mean that files that do not
            // do this are invalid gro files. The failure is with the original writer in that case.
            let pos = [&line[20..28], &line[28..36], &line[36..44]]
                .map(|v| v.trim().parse::<f32>().expect("invalid position float"));

            beads.push(Bead::new(resnum, resname, atomname, atomnum, pos));
        }

        // Finally, read the box vectors.
        let boxvec: [f32; 3] = lines
            .next()
            .expect("unexpected end of file: no box vectors specified")?
            .split_whitespace()
            .take(3)
            .map(|v| {
                v.parse()
                    .expect("invalid box vectors: cannot parse component value")
            })
            .collect::<Vec<_>>()
            .try_into()
            .expect("invalid box vectors: fewer than three components are provided");
        let boxvec = Vec3::from_array(boxvec);

        Ok(Structure { beads, boxvec })
    }

    pub fn write_to_gro(&self, writer: impl Write) -> io::Result<()> {
        let mut writer = BufWriter::new(writer);

        // Write the title.
        writeln!(writer, "solvated by {}", env!("CARGO_PKG_NAME"))?;

        // And the number of atoms in the system.
        let n_atoms = self.beads.len();
        writeln!(writer, "{n_atoms}")?;

        // Write the atoms.
        for chunk in self.beads.chunks(0xffffff) {
            // Format this chunk in parallel.
            let formatted: String = chunk
                .par_iter()
                .fold(
                    || String::new(),
                    |mut s, bead| {
                        bead.format_gro(&mut s);
                        s.push('\n');
                        s
                    },
                )
                .collect();
            // And write it away.
            write!(writer, "{formatted}")?;
        }

        // End with the box vectors.
        let [v1x, v2y, v3z] = self.boxvec.to_array();
        writeln!(writer, "{v1x:.5} {v2y:.5} {v3z:.5}")?;

        writer.flush()
    }

    pub fn read_from_gro_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        Self::read_from_gro(file)
    }

    pub fn write_to_gro_file<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let file = std::fs::File::create(path)?;
        self.write_to_gro(file)
    }
}
