use std::collections::BTreeSet;
use std::io::{self, Read};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use bentopy::core::placement::{Placement, PlacementList};
use bentopy::core::utilities::CLEAR_LINE;
use glam::Vec3;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::args::{Args, Mode, ResnumMode};
use crate::limits::Limits;
use crate::structure::{self, Atom, Molecule, rotate_molecule};

/// Read a placement list to return [`Placements`].
fn read_placement_list(
    placement_list: &str,
    root: Option<PathBuf>,
) -> serde_json::Result<PlacementList> {
    let mut placements: PlacementList = serde_json::from_str(placement_list)?;
    if let Some(root) = root {
        for p in &mut placements.placements {
            if p.path.is_relative() {
                let mut path = root.clone();
                path.push(&p.path);
                p.path = path;
            }
        }
    }
    Ok(placements)
}

impl Mode {
    /// Based on the [`Mode`], prepare the list of molecules.
    ///
    /// If necessary, a molecule is loaded from its structure path, centered, possible its atoms
    /// are filtered down to a subset, or its residues or the whole molecule itself is represented
    /// by single particles.
    ///
    /// The goal is to do as little work as necessary for `Mode`s that aim to make lightweight
    /// representations for inspection.
    fn prepare_molecules(
        &self,
        placements: &[Placement],
        ignore_tags: bool,
        verbose: bool,
    ) -> Result<Vec<Molecule>> {
        let load_molecule = |path| -> Result<_> {
            let prefix = if verbose { "" } else { CLEAR_LINE };
            let suffix = if verbose { "\n" } else { "\r" };
            eprint!("{prefix}\tLoading {path:?}...{suffix}");
            let mut molecule = structure::load_molecule(path)?;
            molecule.translate_to_center();
            Ok(molecule)
        };
        let apply_tag = |molecule: &mut Molecule, tag: Option<&str>| {
            if !ignore_tags && let Some(tag) = tag {
                molecule.apply_tag(tag)
            }
        };

        match self {
            Mode::Full => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    apply_tag(&mut molecule, p.tag());
                    Ok(molecule)
                })
                .collect(),
            Mode::Backbone => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    apply_tag(&mut molecule, p.tag());
                    molecule.atoms.retain(|a| ["N", "CA", "C", "O"].contains(&a.name.as_str()));
                    Ok(molecule)
                })
                .collect(),
            Mode::Alpha => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
                    apply_tag(&mut molecule, p.tag());
                    molecule.atoms.retain(|a| a.name.as_str() == "CA");
                    Ok(molecule)
                })
                .collect(),
            Mode::Residue => placements
                .iter()
                .map(|p| {
                    let mut molecule = load_molecule(&p.path)?;
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
                .collect(),
            Mode::Instance => Ok(placements
                .iter()
                .map(|p| {
                    let mut atom = Atom::dummy();
                    if let Some(tag) = p.tag() {
                        atom.resname = tag.into();
                    }
                    Molecule { atoms: vec![atom] }
                })
                .collect()),
        }
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
    placements: &PlacementList,
    limits: Limits,
    mode: Mode,
    resnum_mode: ResnumMode,
    ignore_tags: bool,
    verbose: bool,
) -> Result<Box<[(String, usize)]>> {
    let min_limits = Vec3::new(
        limits.minx.unwrap_or_default(),
        limits.miny.unwrap_or_default(),
        limits.minz.unwrap_or_default(),
    );
    let size = limits.box_size(Vec3::from(placements.size));
    let placements = &placements.placements;
    // Load the molecules and center them with respect to themselves.
    let molecules = mode.prepare_molecules(placements, ignore_tags, verbose)?;
    writeln!(writer, "{} (v{})", env!("CARGO_BIN_NAME"), bentopy::core::version::VERSION)?;
    let placed = placements.iter().zip(&molecules).map(|(p, m)| {
        let n_placed = p
            .batches
            .iter()
            .map(|b| b.positions.iter().filter(|&&pos| limits.is_inside(pos)).count())
            .sum::<usize>();
        let n_atoms = n_placed * m.atoms.len();
        let name = p.name.clone();
        (name, n_placed, n_atoms)
    });
    let n_atoms: usize = placed.clone().map(|(_, _, n)| n).sum();
    let placed = placed.map(|(name, n_placed, _)| (name, n_placed)).collect();
    let mut instance_resnum = 0;
    writeln!(writer, "{n_atoms}")?;
    for (segidx, (placement, molecule)) in placements.iter().zip(molecules).enumerate() {
        {
            let prefix = if verbose { "" } else { CLEAR_LINE };
            let suffix = if verbose { "\n" } else { "\r" };
            let name = &placement.name;
            let tag = placement.tag().map(|tag| format!(":{tag}")).unwrap_or_default();
            let path = &placement.path;
            eprint!(
                "{prefix}\tNow writing '{name}{tag}' from {path:?} data to output file.{suffix}",
            );
        }
        // Make sure that the printed resnum cannot exceed the 5 allotted columns.
        let segidx = segidx % 100000;

        let mut inserts = Vec::new();
        for atom in &molecule.atoms {
            // HACK: For some reason, converting to string here is necessary for getting correct
            // alignment of ArrayString<U5>.
            let resnum = atom.resnum;
            let resname = atom.resname.to_string();
            let atomname = atom.name.to_string();
            let atomnum = atom.num;
            let insert = match resnum_mode {
                ResnumMode::Instance => format!("{resname:<5}{atomname:>5}{atomnum:>5}"),
                ResnumMode::Segment => format!("{segidx:>5}{resname:<5}{atomname:>5}{atomnum:>5}"),
                ResnumMode::Keep => format!("{resnum:>5}{resname:<5}{atomname:>5}{atomnum:>5}"),
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
            let included_positions: Vec<_> =
                batch.positions.iter().filter(|&&position| limits.is_inside(position)).collect();
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
                            ResnumMode::Keep => {}
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
        if verbose {
            eprintln!("\t\tFormatting took {:.3} s.", t0.elapsed().as_secs_f32());
        }

        let t0 = std::time::Instant::now();
        writer.write_all(output.as_bytes())?;
        writer.flush()?;
        if verbose {
            eprintln!("\t\tWriting took    {:.3} s.", t0.elapsed().as_secs_f32());
        }
    }

    // Write the box vectors.
    let [v1x, v2y, v3z] = size.to_array();
    writeln!(writer, "{v1x:.3} {v2y:.3} {v3z:.3}")?;

    Ok(placed)
}

fn write_cif(
    writer: &mut impl io::Write,
    placements: &PlacementList,
    limits: Limits,
    mode: Mode,
    ignore_tags: bool,
    verbose: bool,
) -> Result<Box<[(String, usize)]>> {
    use std::collections::BTreeMap;

    let placements_slice = &placements.placements;
    let molecules = mode.prepare_molecules(placements_slice, ignore_tags, verbose)?;
    let min_limits = Vec3::new(
        limits.minx.unwrap_or_default(),
        limits.miny.unwrap_or_default(),
        limits.minz.unwrap_or_default(),
    );

    writeln!(writer, "data_bentopy")?;
    writeln!(writer, "#")?;

    let placed: Box<[(String, usize)]> = placements_slice
        .iter()
        .zip(&molecules)
        .map(|(p, _)| {
            let n_placed = p
                .batches
                .iter()
                .map(|b| b.positions.iter().filter(|&&pos| limits.is_inside(pos)).count())
                .sum::<usize>();
            (p.name.clone(), n_placed)
        })
        .collect();
    let entity_asym_ids: Vec<Vec<String>> = molecules
        .iter()
        .enumerate()
        .map(|(idx, molecule)| {
            let entity_id = idx + 1;
            let mut set = BTreeSet::new();
            for atom in &molecule.atoms {
                set.insert(cif_asym_id(entity_id, &atom.chain));
            }
            if set.is_empty() {
                set.insert(cif_asym_id(entity_id, ""));
            }
            set.into_iter().collect()
        })
        .collect();

    writeln!(writer, "loop_")?;
    writeln!(writer, "_entity.id")?;
    writeln!(writer, "_entity.details")?;
    writeln!(writer, "_entity.formula_weight")?;
    writeln!(writer, "_entity.pdbx_description")?;
    writeln!(writer, "_entity.pdbx_ec")?;
    writeln!(writer, "_entity.pdbx_entities_per_biological_unit")?;
    writeln!(writer, "_entity.pdbx_formula_weight_exptl")?;
    writeln!(writer, "_entity.pdbx_formula_weight_exptl_method")?;
    writeln!(writer, "_entity.pdbx_fragment")?;
    writeln!(writer, "_entity.pdbx_modification")?;
    writeln!(writer, "_entity.pdbx_mutation")?;
    writeln!(writer, "_entity.pdbx_number_of_molecules")?;
    writeln!(writer, "_entity.pdbx_parent_entity_id")?;
    writeln!(writer, "_entity.pdbx_target_id")?;
    writeln!(writer, "_entity.src_method")?;
    writeln!(writer, "_entity.type")?;
    for (idx, (placement, (_, n_molecules))) in placements_slice.iter().zip(placed.iter()).enumerate() {
        let entity_id = idx + 1;
        let description = cif_quote(&placement.name);
        let source_filename = placement
            .path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_string();
        let fragment = cif_quote(&source_filename);
        writeln!(
            writer,
            "{entity_id} . 0 {description} . . . . {fragment} . . {n_molecules} . . . polymer"
        )?;
    }
    writeln!(writer, "#")?;

    writeln!(writer, "loop_")?;
    writeln!(writer, "_atom_site.group_PDB")?;
    writeln!(writer, "_atom_site.id")?;
    writeln!(writer, "_atom_site.type_symbol")?;
    writeln!(writer, "_atom_site.label_atom_id")?;
    writeln!(writer, "_atom_site.label_comp_id")?;
    writeln!(writer, "_atom_site.label_asym_id")?;
    writeln!(writer, "_atom_site.label_entity_id")?;
    writeln!(writer, "_atom_site.label_seq_id")?;
    writeln!(writer, "_atom_site.Cartn_x")?;
    writeln!(writer, "_atom_site.Cartn_y")?;
    writeln!(writer, "_atom_site.Cartn_z")?;
    writeln!(writer, "_atom_site.occupancy")?;
    writeln!(writer, "_atom_site.B_iso_or_equiv")?;
    writeln!(writer, "_atom_site.pdbx_PDB_model_num")?;
    let mut atom_id = 1usize;
    for (idx, molecule) in molecules.iter().enumerate() {
        let entity_id = idx + 1;
        for atom in &molecule.atoms {
            let asym_id = cif_asym_id(entity_id, &atom.chain);
            let x = atom.pos.x as f64 * 10.0;
            let y = atom.pos.y as f64 * 10.0;
            let z = atom.pos.z as f64 * 10.0;
            writeln!(
                writer,
                "ATOM {atom_id} C {} {} {asym_id} {entity_id} {} {x:.3} {y:.3} {z:.3} 1.00 0.00 1",
                atom.name,
                atom.resname,
                atom.resnum
            )?;
            atom_id += 1;
        }
    }
    writeln!(writer, "#")?;

    #[derive(Clone)]
    struct Operator {
        id: usize,
        rot: [[f32; 3]; 3],
        tr: [f32; 3],
    }

    let mut op_ids_by_transform: BTreeMap<String, usize> = BTreeMap::new();
    let mut operators: Vec<Operator> = Vec::new();
    let mut per_entity_operator_ids: Vec<Vec<usize>> = Vec::new();
    for (placement, molecule) in placements_slice.iter().zip(&molecules) {
        let mut entity_ops = Vec::new();
        for batch in &placement.batches {
            let rot = batch.rotation;
            let rotated = rotate_molecule(molecule, &batch.rotation());
            let min = rotated.min() + min_limits;
            let m = row_major_to_mmcif_matrix(rot);
            for tr in &batch.positions {
                if !limits.is_inside(*tr) {
                    continue;
                }
                let tr = [tr[0] - min.x, tr[1] - min.y, tr[2] - min.z];
                let key = format!(
                    "{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}|{:.6},{:.6},{:.6}",
                    m[0],
                    m[1],
                    m[2],
                    m[3],
                    m[4],
                    m[5],
                    m[6],
                    m[7],
                    m[8],
                    tr[0],
                    tr[1],
                    tr[2]
                );
                let id = if let Some(id) = op_ids_by_transform.get(&key) {
                    *id
                } else {
                    let id = op_ids_by_transform.len() + 1;
                    op_ids_by_transform.insert(key, id);
                    operators.push(Operator { id, rot, tr });
                    id
                };
                entity_ops.push(id);
            }
        }
        entity_ops.sort_unstable();
        entity_ops.dedup();
        per_entity_operator_ids.push(entity_ops);
    }
    operators.sort_by_key(|op| op.id);

    writeln!(writer, "loop_")?;
    writeln!(writer, "_pdbx_struct_assembly.id")?;
    writeln!(writer, "_pdbx_struct_assembly.details")?;
    writeln!(writer, "_pdbx_struct_assembly.method_details")?;
    writeln!(writer, "_pdbx_struct_assembly.oligomeric_details")?;
    writeln!(writer, "_pdbx_struct_assembly.oligomeric_count")?;
    let n_instances: usize = placed.iter().map(|(_, n)| *n).sum();
    writeln!(writer, "1 'author_and_software_defined_assembly' bentopy ? {n_instances}")?;
    writeln!(writer, "#")?;

    writeln!(writer, "loop_")?;
    writeln!(writer, "_pdbx_struct_assembly_gen.assembly_id")?;
    writeln!(writer, "_pdbx_struct_assembly_gen.oper_expression")?;
    writeln!(writer, "_pdbx_struct_assembly_gen.asym_id_list")?;
    for (idx, operator_ids) in per_entity_operator_ids.iter().enumerate() {
        let asym_id_list = entity_asym_ids[idx].join(",");
        let oper_expression = cif_quote(&format!("({})", compress_ranges(operator_ids)));
        writeln!(writer, "1 {oper_expression} {asym_id_list}")?;
    }
    writeln!(writer, "#")?;

    writeln!(writer, "loop_")?;
    writeln!(writer, "_pdbx_struct_oper_list.id")?;
    writeln!(writer, "_pdbx_struct_oper_list.type")?;
    writeln!(writer, "_pdbx_struct_oper_list.name")?;
    writeln!(writer, "_pdbx_struct_oper_list.symmetry_operation")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[1][1]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[1][2]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[1][3]")?;
    writeln!(writer, "_pdbx_struct_oper_list.vector[1]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[2][1]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[2][2]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[2][3]")?;
    writeln!(writer, "_pdbx_struct_oper_list.vector[2]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[3][1]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[3][2]")?;
    writeln!(writer, "_pdbx_struct_oper_list.matrix[3][3]")?;
    writeln!(writer, "_pdbx_struct_oper_list.vector[3]")?;
    for op in operators {
        let op_type = if op.id == 1 { "identity operation" } else { "point symmetry operation" };
        let m = row_major_to_mmcif_matrix(op.rot);
        writeln!(
            writer,
            "{id} '{}' ? ? {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
            op_type,
            m[0],
            m[1],
            m[2],
            op.tr[0] as f64 * 10.0,
            m[3],
            m[4],
            m[5],
            op.tr[1] as f64 * 10.0,
            m[6],
            m[7],
            m[8],
            op.tr[2] as f64 * 10.0,
            id = op.id
        )?;
    }
    writeln!(writer, "#")?;

    Ok(placed)
}

fn cif_quote(value: &str) -> String {
    let escaped = value.replace('\'', "''");
    format!("'{escaped}'")
}

fn cif_asym_id(idx: usize, chain: &str) -> String {
    let base = cif_entity_code(idx);
    let chain = chain.trim();
    if chain.is_empty() {
        format!("{base}A")
    } else {
        format!("{base}{chain}")
    }
}

fn cif_entity_code(mut idx: usize) -> String {
    const ALPH: &[u8; 26] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let mut out = Vec::new();
    while idx > 0 {
        let rem = (idx - 1) % 26;
        out.push(ALPH[rem] as char);
        idx = (idx - 1) / 26;
    }
    out.reverse();
    out.into_iter().collect()
}

fn compress_ranges(ids: &[usize]) -> String {
    if ids.is_empty() {
        return "1".to_string();
    }

    let mut parts = Vec::new();
    let mut start = ids[0];
    let mut prev = ids[0];
    for &id in &ids[1..] {
        if id == prev + 1 {
            prev = id;
            continue;
        }
        if start == prev {
            parts.push(start.to_string());
        } else {
            parts.push(format!("{start}-{prev}"));
        }
        start = id;
        prev = id;
    }

    if start == prev {
        parts.push(start.to_string());
    } else {
        parts.push(format!("{start}-{prev}"));
    }

    parts.join(",")
}

fn row_major_to_mmcif_matrix(rot: [[f32; 3]; 3]) -> [f32; 9] {
    // bentopy stores batch.rotation as RowMajorRotation.
    // mmCIF _pdbx_struct_oper_list.matrix[i][j] is also row/column indexed, so no transpose.
    [
        rot[0][0],
        rot[0][1],
        rot[0][2],
        rot[1][0],
        rot[1][1],
        rot[1][2],
        rot[2][0],
        rot[2][1],
        rot[2][2],
    ]
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
        env!("CARGO_BIN_NAME"),
        bentopy::core::version::VERSION
    )?;
    writeln!(writer, "; {}", env!("CARGO_PKG_REPOSITORY"))?;
    writeln!(writer, ";")?;
    let citation = bentopy::core::citation::CITATION.lines();
    for line in citation {
        writeln!(writer, "; {line}")?;
    }
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
pub fn render(args: Args) -> anyhow::Result<()> {
    let Args {
        input: input_path,
        output: output_path,
        topol: topol_path,
        root,
        limits,
        mode,
        resnum_mode,
        ignore_tags,
        verbose,
    } = args;

    eprint!("Reading from {input_path:?}... ");
    let t0 = std::time::Instant::now();
    let mut placement_list = String::new();
    if input_path.to_str() == Some("-") {
        std::io::stdin().read_to_string(&mut placement_list)?;
    } else {
        std::fs::File::open(&input_path)?.read_to_string(&mut placement_list)?;
    }
    let placements = read_placement_list(&placement_list, root)
        .with_context(|| format!("Failed to read placement list from {input_path:?}"))?;
    eprintln!("Done in {:.3} ms.", t0.elapsed().as_millis());

    eprintln!("Writing structure to {output_path:?}... ");
    let t0 = std::time::Instant::now();
    let outfile = std::fs::File::create(&output_path)?;
    let limits = limits.unwrap_or_default();
    let extension = output_path.extension().and_then(|s| s.to_str()).map(|s| s.to_ascii_lowercase());
    let placed = match extension.as_deref() {
        Some("gro") => write_gro(
            &mut std::io::BufWriter::new(outfile),
            &placements,
            limits,
            mode,
            resnum_mode,
            ignore_tags,
            verbose,
        )?,
        Some("cif") | Some("mmcif") => write_cif(
            &mut std::io::BufWriter::new(outfile),
            &placements,
            limits,
            mode,
            ignore_tags,
            verbose,
        )?,
        _ => bail!("Unsupported output extension for {output_path:?}. Use .gro, .cif, or .mmcif."),
    };
    eprintln!("{CLEAR_LINE}Done in {:.3} s.", t0.elapsed().as_secs_f32());
    eprintln!(
        "Wrote {} placements of {} different kinds of structures.",
        placed.iter().map(|(_, c)| c).sum::<usize>(),
        placed.iter().filter(|(_, c)| *c > 0).count(),
    );

    if let Some(topol_path) = topol_path {
        let mut top = std::fs::File::create(&topol_path)?;
        write_top(&mut top, &placements.topol_includes, &placements.title, &placed)?;
        eprintln!("Wrote topology to {topol_path:?}.");
    }

    Ok(())
}
