use std::io::{self, Read};
use std::path::PathBuf;
use std::str::FromStr;

use clap::ValueEnum;
use glam::{Mat3, Vec3};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::Deserialize;

use crate::limits::Limits;
use crate::structure::{load_molecule, rotate_molecule, Atom, Molecule};

/// A space with a set of [`Placement`]s.
#[derive(Debug, Clone, Deserialize)]
struct Placements {
    title: String,
    size: [f32; 3],
    topol_includes: Vec<String>,
    placements: Vec<Placement>,
}

/// A set of positions associated with some structure.
///
/// The structure is named and is stored at the provided `path`.
///
/// Positions are stored in [`Batch`]es.
#[derive(Debug, Clone, Deserialize)]
struct Placement {
    name: String,
    /// A tag that will replace the associated structure's `resname` field, if present.
    tag: Option<String>,
    /// Path to the segment's molecule file.
    path: PathBuf,
    batches: Vec<Batch>,
}

impl Placement {
    fn tag(&self) -> Option<&str> {
        self.tag.as_ref().map(|t| t.as_str())
    }
}

/// A set of positions that share a specific rotation.
#[derive(Debug, Clone, Deserialize)]
struct Batch {
    /// The 3×3 rotation matrix.
    ///
    /// The data is stored in a row-major order.
    rotation: [[f32; 3]; 3],
    positions: Vec<[f32; 3]>,
}

impl Batch {
    /// Returns the 3×3 rotation matrix of this [`Batch`].
    fn rotation(&self) -> Mat3 {
        // Because the batch rotation is stored in a row-major order, and the initializer we use
        // here assumes column-major order, we need to transpose the matrix before using it.
        Mat3::from_cols_array_2d(&self.rotation).transpose()
    }
}

/// Read a placement list to return [`Placements`].
fn read_placement_list(placement_list: &str, root: Option<PathBuf>) -> Placements {
    let mut placements: Placements = serde_json::from_str(placement_list).unwrap();
    if let Some(root) = root {
        for p in &mut placements.placements {
            if p.path.is_relative() {
                let mut path = root.clone();
                path.push(&p.path);
                p.path = path;
            }
        }
    }
    placements
}

#[derive(Debug, Default, Clone, ValueEnum)]
pub enum ResnumMode {
    /// Each segment instance has a unique residue number.
    #[default]
    Instance,
    /// All instances of a segment have the same residue number.
    Segment,
}

impl FromStr for ResnumMode {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "instance" => Self::Instance,
            "segment" => Self::Segment,
            _ => return Err(()),
        })
    }
}

/// Write out [`Placements`] as a `gro` file.
///
/// See the Gromacs manual entry on the [gro file][gro_manual].
///
/// Note that each [`Placement`] (i.e., kind of structure) is given a different `resnum`, such that
/// the different kinds of structures can be targeted individually in programmes that use this
/// `gro` file output. For example, each kind of molecule may be drawn in a different color in a
/// molecular visualization program.
///
/// [gro_manual]: https://manual.gromacs.org/archive/5.0.3/online/gro.html
fn write_gro(
    writer: &mut impl io::Write,
    placements: &Placements,
    limits: Limits,
    mode: Mode,
    resnum_mode: ResnumMode,
    ignore_tags: bool,
) -> io::Result<Box<[(String, usize)]>> {
    let apply_tag = |molecule: &mut Molecule, tag: Option<&str>| {
        if !ignore_tags {
            if let Some(tag) = tag {
                molecule.apply_tag(tag)
            }
        }
    };

    let min_limits = Vec3::new(
        limits.minx.unwrap_or_default(),
        limits.miny.unwrap_or_default(),
        limits.minz.unwrap_or_default(),
    );
    let size = limits.box_size(Vec3::from(placements.size));
    let placements = &placements.placements;
    // Load the molecules and center them with respect to themselves.
    let molecules: Vec<_> = {
        match mode {
            Mode::Full => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    apply_tag(&mut molecule, p.tag());
                    molecule.translate_to_center();
                    Ok(molecule)
                })
                .collect::<io::Result<Vec<_>>>()?,
            Mode::Backbone => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    apply_tag(&mut molecule, p.tag());
                    molecule.translate_to_center();
                    molecule.atoms = molecule
                        .atoms
                        .into_iter()
                        .filter(|a| ["N", "CA", "C", "O"].contains(&a.name.as_str()))
                        .collect();
                    Ok(molecule)
                })
                .collect::<io::Result<Vec<_>>>()?,
            Mode::Alpha => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    apply_tag(&mut molecule, p.tag());
                    molecule.translate_to_center();
                    molecule.atoms = molecule
                        .atoms
                        .into_iter()
                        .filter(|a| a.name.as_str() == "CA")
                        .collect();
                    Ok(molecule)
                })
                .collect::<io::Result<Vec<_>>>()?,
            Mode::Residue => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    molecule.translate_to_center();
                    let mut residues = Vec::new();
                    let mut residue_atoms = Vec::new();
                    for atom in molecule.atoms {
                        let Some(last) = residue_atoms.last() else {
                            residue_atoms.push(atom);
                            continue;
                        };
                        if last.resnum == atom.resnum {
                            residue_atoms.push(atom);
                        } else {
                            // The residue we just collected is complete.
                            let avg_pos =
                                residue_atoms.iter().fold(Vec3::ZERO, |acc, a| acc + a.pos)
                                    / residue_atoms.len() as f32;
                            residues.push(Atom {
                                name: "DUMMY".into(),
                                num: residues.len() as u32,
                                pos: avg_pos,
                                ..last.clone()
                            });

                            // Prepare for the next resdue.
                            residue_atoms.clear();
                            residue_atoms.push(atom);
                        }
                    }
                    molecule.atoms = residues;
                    apply_tag(&mut molecule, p.tag());
                    molecule.translate_to_center();
                    Ok(molecule)
                })
                .collect::<io::Result<Vec<_>>>()?,
            Mode::Instance => placements
                .iter()
                .map(|p| {
                    let mut atom = Atom::dummy();
                    if let Some(tag) = p.tag() {
                        atom.resname = tag.into();
                    }
                    Molecule { atoms: vec![atom] }
                })
                .collect(),
        }
    };
    writeln!(writer, "{}", env!("CARGO_PKG_NAME"))?;
    let placed = placements.iter().zip(&molecules).map(|(p, m)| {
        let n_placed = p
            .batches
            .iter()
            .map(|b| {
                b.positions
                    .iter()
                    .filter(|&&pos| limits.is_inside(pos))
                    .count()
            })
            .sum::<usize>();
        let n_atoms = n_placed * m.atoms.len();
        let name = p.name.clone();
        (name, n_placed, n_atoms)
    });
    let n_atoms: usize = placed.clone().map(|(_, _, n)| n).sum();
    let placed = placed.map(|(name, n_placed, _)| (name, n_placed)).collect();
    let mut instance_resnum = 0;
    writeln!(writer, "{n_atoms}")?;
    for (resnum, (placement, molecule)) in placements.iter().zip(molecules).enumerate() {
        eprintln!(
            "\tNow writing '{}' from {:?} data to output file.",
            placement.name, placement.path
        );
        // Make sure that the printed resnum cannot exceed the 5 allotted columns.
        let resnum = resnum % 100000;

        let mut inserts = Vec::new();
        for atom in &molecule.atoms {
            // HACK: For some reason, converting to string here is necessary for getting correct
            // alignment of ArrayString<U5>.
            let resname = atom.resname.to_string();
            let atomname = atom.name.to_string();
            let atomnum = atom.num;
            let insert = match resnum_mode {
                ResnumMode::Instance => format!("{resname:<5}{atomname:>5}{atomnum:>5}"),
                ResnumMode::Segment => format!("{resnum:>5}{resname:<5}{atomname:>5}{atomnum:>5}"),
            };
            inserts.push(insert)
        }

        let mut output = String::new();
        let t0 = std::time::Instant::now();
        for batch in &placement.batches {
            // This tends to take around half a ms.
            let molecule = rotate_molecule(&molecule, &batch.rotation());

            // Move over the positions so that they do not float somewhere in space, far from the
            // origin, in case the limits are set. This puts them within the specified limits.
            let min = molecule.min() + min_limits;
            let included_positions: Vec<_> = batch
                .positions
                .iter()
                .filter(|&&position| limits.is_inside(position))
                .collect();
            let n_included_instances = included_positions.len();
            let out: String = included_positions
                .par_iter()
                .enumerate()
                .map(|(batch_instance_resnum, position)| {
                    let mut output = String::new();
                    for (atom, insert) in molecule.atoms.iter().zip(&inserts) {
                        let mut pos = atom.pos;
                        pos.x += position[0] - min.x;
                        pos.y += position[1] - min.y;
                        pos.z += position[2] - min.z;

                        match resnum_mode {
                            ResnumMode::Instance => {
                                let resnum = instance_resnum + batch_instance_resnum;
                                // Make sure that the printed resnum cannot exceed the 5 allotted columns.
                                let resnum = resnum % 100000;
                                output.push_str(&format!("{resnum:>5}"))
                            }
                            ResnumMode::Segment => {}
                        };
                        output.push_str(insert);
                        output.push_str(&format!("{:8.3}{:8.3}{:8.3}\n", pos.x, pos.y, pos.z));
                    }
                    output
                })
                .collect();
            instance_resnum += n_included_instances;
            output.push_str(&out);
        }
        let dt = std::time::Instant::now() - t0;
        eprintln!("\t\tFormatting took {:.3} s.", dt.as_secs_f32());

        let t0 = std::time::Instant::now();
        writer.write_all(output.as_bytes())?;
        writer.flush()?;
        let dt = std::time::Instant::now() - t0;
        eprintln!("\t\tWriting took    {:.3} s.", dt.as_secs_f32());
    }

    // Write the box vectors.
    let [v1x, v2y, v3z] = size.to_array();
    writeln!(writer, "{v1x:.3} {v2y:.3} {v3z:.3}")?;

    Ok(placed)
}

#[derive(Debug, Default, Clone, ValueEnum)]
pub enum Mode {
    /// Output all atoms.
    #[default]
    Full,
    /// Only the backbone atoms.
    Backbone,
    /// Only the alpha carbons.
    Alpha,
    /// One bead per residue.
    Residue,
    /// One bead per instance of a structure.
    Instance,
}

impl FromStr for Mode {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "full" => Self::Full,
            "backbone" | "bb" => Self::Full,
            "alpha" | "ca" | "a" => Self::Alpha,
            "residue" | "res" => Self::Residue,
            "instance" => Self::Instance,
            _ => return Err(()),
        })
    }
}

/// Write out the topology as a `top` file.
///
/// See the Gromacs manual entry on the [top file][top_manual].
///
/// [top_manual]: https://manual.gromacs.org/archive/5.0/online/top.html
fn write_top(
    writer: &mut impl io::Write,
    includes: &[String],
    title: &str,
    placed: &[(String, usize)],
) -> io::Result<()> {
    // Write some meta stuff.
    writeln!(
        writer,
        "; This topology file was generated by {} (v{}).",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(writer, "; {}, 2024.", env!("CARGO_PKG_AUTHORS"))?;
    writeln!(writer, "; {}", env!("CARGO_PKG_REPOSITORY"))?;
    writeln!(writer)?;

    // Dump the includes.
    for include in includes {
        writeln!(writer, "#include \"{include}\"")?;
    }
    writeln!(writer)?;

    // Write the system title.
    writeln!(writer, "[ system ]")?;
    writeln!(writer, "{title}")?;
    writeln!(writer)?;

    // Write out the name of each structure and how many of them were written to the output file.
    writeln!(writer, "[ molecules ]")?;
    for (name, count) in placed {
        if *count > 0 {
            writeln!(writer, "{name}\t{count}")?;
        }
    }

    Ok(())
}

/// Render a placement list to a gro file.
///
/// If `input_path` is "-", the placement list is read from stdin.
pub fn render(
    input_path: PathBuf,
    output_path: PathBuf,
    topol_path: Option<PathBuf>,
    root: Option<PathBuf>,
    limits: Option<Limits>,
    mode: Mode,
    resnum_mode: ResnumMode,
    ignore_tags: bool,
) -> io::Result<()> {
    eprint!("Reading from {input_path:?}... ");
    let t0 = std::time::Instant::now();
    let mut placement_list = String::new();
    if input_path.to_str() == Some("-") {
        std::io::stdin().read_to_string(&mut placement_list)?;
    } else {
        std::fs::File::open(&input_path)?.read_to_string(&mut placement_list)?;
    }
    let placements = read_placement_list(&placement_list, root);
    let dt = std::time::Instant::now() - t0;
    eprintln!("Done in {:.3} ms.", dt.as_millis());

    eprintln!("Writing structure to {output_path:?}... ");
    let t0 = std::time::Instant::now();
    let outfile = std::fs::File::create(&output_path)?;
    let placed = write_gro(
        &mut std::io::BufWriter::new(outfile),
        &placements,
        limits.unwrap_or_default(),
        mode,
        resnum_mode,
        ignore_tags,
    )?;
    let dt = std::time::Instant::now() - t0;
    eprintln!("Done in {:.3} s.", dt.as_secs_f32());
    eprintln!(
        "Wrote {} placements of {} different kinds of structures.",
        placed.iter().map(|(_, c)| c).sum::<usize>(),
        placed.iter().filter(|(_, c)| *c > 0).count(),
    );

    if let Some(topol_path) = topol_path {
        let mut top = std::fs::File::create(&topol_path)?;
        write_top(
            &mut top,
            &placements.topol_includes,
            &placements.title,
            &placed,
        )?;
        eprintln!("Wrote topology to {topol_path:?}.");
    }

    Ok(())
}
