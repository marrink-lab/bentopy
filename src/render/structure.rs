use std::{collections::HashMap, io::Read, path::Path};

use anyhow::{Context, Result};
use eightyseven::structure::{AtomName, AtomNum, ResName, ResNum};
use glam::{Mat3, Vec3};

/// A molecule as derived from PDB 'ATOM' records.
#[derive(Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl Molecule {
    /// Create a new molecule from a PDB file.
    fn from_pdb(pdb: &str) -> Result<Self> {
        let mut atoms = Vec::new();
        for (ln, line) in pdb.lines().enumerate() {
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let ln = ln + 1;
                let atom = Atom::from_pdb_atom_line(line)
                    .with_context(|| format!("Could not parse atom on line {ln}"))?;
                atoms.push(atom);
            }
        }

        Ok(Self { atoms })
    }

    /// Create a new molecule from a gro file.
    fn from_gro(gro: &str) -> Result<Self> {
        let mut atoms = Vec::new();
        let mut lines = gro.lines().enumerate();
        let (_ln, _title) = lines.next().context("Expected a title line")?;
        let n_atoms: usize = lines
            .next()
            .context("Expected the number of atoms")?
            .1
            .trim()
            .parse()
            .context("Could not parse the number of atoms")?;
        // Read the atoms.
        for (ln, line) in lines.take(n_atoms) {
            let ln = ln + 1;
            let atom = Atom::from_gro_atom_line(line)
                .with_context(|| format!("Could not parse atom on line {ln}"))?;
            atoms.push(atom);
        }
        // We don't check for the presence and correctness of the box vectors, even though they
        // should be there. We just don't care here.

        Ok(Self { atoms })
    }

    /// Create a new molecule from an mmCIF file.
    fn from_cif(cif: &str) -> Result<Self> {
        let parsed = parse_cif_data(cif)?;
        if parsed.operations.is_empty() || parsed.assembly_rules.is_empty() {
            return Ok(Self { atoms: parsed.flat_atoms() });
        }
        Ok(Self { atoms: build_cif_assembly(&parsed, "1")? })
    }

    /// Apply a `tag` to the atoms in this [`Molecule`].
    ///
    /// This will replace the `resname` for each atom with the provided `tag`.
    pub fn apply_tag(&mut self, tag: impl Into<ResName>) {
        let tag = tag.into();
        for Atom { resname, .. } in &mut self.atoms {
            *resname = tag;
        }
    }

    /// Returns the center of this [`Molecule`].
    ///
    /// Here, center refers to mean position. This can be understood as a center of mass where all
    /// atoms are considered to have the same weight.
    pub fn center(&self) -> Vec3 {
        if self.atoms.is_empty() {
            return Vec3::ZERO;
        }
        let mut mean = Vec3::ZERO;
        for atom in &self.atoms {
            mean += atom.pos
        }
        mean / self.atoms.len() as f32
    }

    /// Translate all atoms such that their center lies at the origin.
    ///
    /// Here, center refers to mean position. This can be understood as a center of mass where all
    /// atoms are considered to have the same weight.
    pub fn translate_to_center(&mut self) {
        let center = self.center();
        for atom in &mut self.atoms {
            atom.pos -= center;
        }
    }

    /// Returns the minimum _x_, _y_, and _z_ values of this [`Molecule`] as a [`Vec3`].
    ///
    /// Note that this is a virtual point, not an actual point in the atoms point cloud.
    /// In other words, it is not the lowest point in a sorted list of existing positions.
    pub fn min(&self) -> Vec3 {
        let mut min = Vec3::ZERO;
        for atom in &self.atoms {
            let pos = atom.pos;
            if pos.x < min.x {
                min.x = pos.x
            }
            if pos.y < min.y {
                min.y = pos.y
            }
            if pos.z < min.z {
                min.z = pos.z
            }
        }
        min
    }
}

/// Assumes the [Molecule]'s atoms have been centered.
pub fn rotate_molecule(molecule: &Molecule, rotmat: &Mat3) -> Molecule {
    let mut molecule = molecule.clone();
    for atom in &mut molecule.atoms {
        atom.pos = rotmat.mul_vec3(atom.pos)
    }

    molecule
}

/// The representation of an atom.
#[derive(Clone)]
pub struct Atom {
    /// Atom name.
    pub name: AtomName,
    /// Residue name.
    pub resname: ResName,
    /// Residue sequence number.
    pub resnum: ResNum,
    /// Number of the atom within its structure.
    ///
    /// Also understood as its 'serial' number.
    pub num: AtomNum,
    /// Position in nanometers.
    pub pos: Vec3,
    /// Chain identifier.
    pub chain: String,
}

impl Atom {
    // COLUMNS        DATA  TYPE    FIELD        DEFINITION
    // -------------------------------------------------------------------------------------
    //  1 -  6        Record name   "ATOM  "
    //  7 - 11        Integer       serial       Atom  serial number.
    // 13 - 16        Atom          name         Atom name.
    // 17             Character     altLoc       Alternate location indicator.
    // 18 - 20        Residue name  resName      Residue name.
    // 22             Character     chainID      Chain identifier.
    // 23 - 26        Integer       resSeq       Residue sequence number.
    // 27             AChar         iCode        Code for insertion of residues.
    // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    // 55 - 60        Real(6.2)     occupancy    Occupancy.
    // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    // 77 - 78        LString(2)    element      Element symbol, right-justified.
    // 79 - 80        LString(2)    charge       Charge  on the atom.
    //
    // Examples:
    // ATOM      1  N   ALA A   1      11.104   6.134  -6.504
    // ATOM      1  N   MET     1      48.048  69.220  58.803
    // ATOM      2  C13 DPPCM 367      31.671 -46.874  39.426  1.00  0.00      MEMB
    /// Read a single "ATOM" or "HETATM" record from a PDB and return an Atom.
    fn from_pdb_atom_line(line: &str) -> Result<Atom> {
        let serial = line[6..11].trim().parse().context("Could not parse atom serial")?;
        let name = line[12..16].trim().into();
        // NOTE: Even though the PDB specification only regards columns 18..21 as constituting the
        // resname, in practice the character directly after that is also included. This column is
        // not defined by the spec. Especially for telling apart lipids like DPPC and DPPG, it's
        // quite important to include that by-convention resname character. Thanks Abby!
        let resname = line[17..21].trim().into();
        let resnum = line[22..26].trim().parse().context("Could not parse atom resnum")?;
        let x = line[30..38].trim().parse().context("Could not parse x coordinate")?;
        let y = line[38..46].trim().parse().context("Could not parse y coordinate")?;
        let z = line[46..54].trim().parse().context("Could not parse z coordinate")?;
        Ok(Atom {
            name,
            resname,
            resnum,
            num: serial,
            pos: Vec3::new(x, y, z) / 10.0, // Convert from Å to nm.
            chain: line[21..22].trim().to_string(),
        })
    }

    // `    2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244`
    /// Read a single atom line from a gro file and return an Atom.
    fn from_gro_atom_line(line: &str) -> Result<Atom> {
        let resnum = line[0..5].trim().parse().context("Could not parse resnum")?;
        let resname = line[5..10].trim().into();
        let name = line[10..15].trim().into(); // Atom name.
        let num = line[15..20].trim().parse().context("Could not parse atomnum")?; // Atom number.
        let x = line[20..28].trim().parse().context("Could not parse x coordinate")?;
        let y = line[28..36].trim().parse().context("Could not parse y coordinate")?;
        let z = line[36..44].trim().parse().context("Could not parse z coordinate")?;
        Ok(Atom {
            name,
            resname,
            resnum,
            num,
            pos: Vec3::new(x, y, z), // Values are already in nm.
            chain: String::new(),
        })
    }

    fn from_cif_atom_fields(
        fields: &[&str],
        columns: &CifAtomSiteColumns,
        fallback_serial: u32,
    ) -> Result<Atom> {
        let x_idx = columns.x_idx.context("Missing _atom_site.Cartn_x in mmCIF loop")?;
        let y_idx = columns.y_idx.context("Missing _atom_site.Cartn_y in mmCIF loop")?;
        let z_idx = columns.z_idx.context("Missing _atom_site.Cartn_z in mmCIF loop")?;

        let x: f32 = fields.get(x_idx).context("Missing x coordinate")?.parse().context("Could not parse x coordinate")?;
        let y: f32 = fields.get(y_idx).context("Missing y coordinate")?.parse().context("Could not parse y coordinate")?;
        let z: f32 = fields.get(z_idx).context("Missing z coordinate")?.parse().context("Could not parse z coordinate")?;

        let num: AtomNum = match columns.id_idx.and_then(|idx| fields.get(idx).copied()) {
            Some(id) if id != "." && id != "?" => id.parse().unwrap_or(fallback_serial),
            _ => fallback_serial,
        };

        let name: AtomName = columns
            .name_idx
            .and_then(|idx| fields.get(idx).copied())
            .filter(|v| *v != "." && *v != "?")
            .unwrap_or("DUMMY")
            .into();

        let resname: ResName = columns
            .resname_idx
            .and_then(|idx| fields.get(idx).copied())
            .filter(|v| *v != "." && *v != "?")
            .unwrap_or("DUM")
            .into();

        let resnum: ResNum = match columns.resnum_idx.and_then(|idx| fields.get(idx).copied()) {
            Some(v) if v != "." && v != "?" => v.parse().unwrap_or(0),
            _ => 0,
        };
        let chain = columns
            .chain_idx
            .and_then(|idx| fields.get(idx).copied())
            .filter(|v| *v != "." && *v != "?")
            .unwrap_or("")
            .to_string();

        Ok(Atom {
            name,
            resname,
            resnum,
            num,
            pos: Vec3::new(x, y, z) / 10.0, // Convert from A to nm.
            chain,
        })
    }

    /// Create a dummy atom.
    pub fn dummy() -> Self {
        Self {
            name: "DUMMY".into(),
            resname: "DUMMY".into(),
            resnum: 0,
            num: 0,
            pos: Vec3::ZERO,
            chain: String::new(),
        }
    }
}

/// Load a [`Molecule`] from a pdb file.
pub fn load_molecule<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Molecule> {
    let mut data = String::new();
    std::fs::File::open(&path)?.read_to_string(&mut data)?;

    let molecule = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gro") => Molecule::from_gro(&data)
            .with_context(|| format!("Could not parse .gro file at {path:?}"))?,
        Some("pdb") => Molecule::from_pdb(&data)
            .with_context(|| format!("Could not parse .pdb file at {path:?}"))?,
        Some("cif") | Some("mmcif") => Molecule::from_cif(&data)
            .with_context(|| format!("Could not parse .cif file at {path:?}"))?,
        None | Some(_) => {
            eprintln!("WARNING: Assuming {path:?} is a pdb file.");
            Molecule::from_pdb(&data)
                .with_context(|| format!("Could not parse file at {path:?} as .pdb"))?
        }
    };

    Ok(molecule)
}

fn split_cif_tokens(line: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let bytes = line.as_bytes();
    let mut i = 0;

    while i < bytes.len() {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }

        let start = i;
        let quote = if bytes[i] == b'\'' || bytes[i] == b'"' {
            let q = bytes[i];
            i += 1;
            while i < bytes.len() && bytes[i] != q {
                i += 1;
            }
            if i < bytes.len() {
                i += 1;
            }
            Some(q)
        } else {
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            None
        };

        let token = &line[start..i];
        if quote.is_some() && token.len() >= 2 {
            tokens.push(&token[1..token.len() - 1]);
        } else {
            tokens.push(token);
        }
    }

    tokens
}

#[derive(Clone, Copy)]
struct Transform {
    rot: [[f32; 3]; 3],
    tr: Vec3,
}

impl Transform {
    fn identity() -> Self {
        Self { rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], tr: Vec3::ZERO }
    }
}

#[derive(Clone)]
struct AssemblyRule {
    assembly_id: String,
    asym_ids: Vec<String>,
    oper_expression: String,
}

#[derive(Default)]
struct ParsedCifData {
    atoms_by_asym: HashMap<String, Vec<Atom>>,
    operations: HashMap<String, Transform>,
    assembly_rules: Vec<AssemblyRule>,
}

impl ParsedCifData {
    fn flat_atoms(&self) -> Vec<Atom> {
        let mut out = Vec::new();
        for atoms in self.atoms_by_asym.values() {
            out.extend_from_slice(atoms);
        }
        out
    }
}

fn parse_cif_data(cif: &str) -> Result<ParsedCifData> {
    let mut parsed = ParsedCifData::default();
    let lines: Vec<(usize, &str)> = cif.lines().enumerate().map(|(i, l)| (i + 1, l)).collect();
    let mut i = 0usize;
    while i < lines.len() {
        let t = lines[i].1.trim();
        if t != "loop_" {
            i += 1;
            continue;
        }
        let (headers, rows, next_i) = parse_cif_loop_rows(&lines, i + 1);
        i = next_i;
        if headers.is_empty() {
            continue;
        }
        let first = headers[0].as_str();
        if first.starts_with("_atom_site.") {
            parse_atom_site_loop_rows(&headers, &rows, &mut parsed)?;
        } else if first.starts_with("_pdbx_struct_oper_list.") {
            parse_oper_list_loop_rows(&headers, &rows, &mut parsed)?;
        } else if first.starts_with("_pdbx_struct_assembly_gen.") {
            parse_assembly_gen_loop_rows(&headers, &rows, &mut parsed)?;
        }
    }
    Ok(parsed)
}

fn parse_cif_loop_rows(
    lines: &[(usize, &str)],
    mut i: usize,
) -> (Vec<String>, Vec<(usize, Vec<String>)>, usize) {
    let mut headers = Vec::new();
    while i < lines.len() {
        let t = lines[i].1.trim();
        if t.starts_with('_') {
            if let Some(h) = t.split_whitespace().next() {
                headers.push(h.to_string());
            }
            i += 1;
        } else {
            break;
        }
    }

    let mut rows = Vec::new();
    let ncols = headers.len();
    let mut buffer: Vec<String> = Vec::new();
    let mut row_start_ln: Option<usize> = None;
    while i < lines.len() {
        let (ln, raw) = lines[i];
        let t = raw.trim();
        if t.is_empty() || t == "#" {
            i += 1;
            continue;
        }
        if t == "loop_" || t == "stop_" || t.starts_with("data_") || t.starts_with('_') {
            break;
        }
        let tokens: Vec<String> = split_cif_tokens(t).into_iter().map(str::to_string).collect();
        if !tokens.is_empty() {
            if row_start_ln.is_none() {
                row_start_ln = Some(ln);
            }
            buffer.extend(tokens);
            while ncols > 0 && buffer.len() >= ncols {
                let fields: Vec<String> = buffer.drain(0..ncols).collect();
                rows.push((row_start_ln.unwrap_or(ln), fields));
                row_start_ln = if buffer.is_empty() { None } else { Some(ln) };
            }
        }
        i += 1;
    }
    (headers, rows, i)
}

fn parse_atom_site_loop_rows(
    headers: &[String],
    rows: &[(usize, Vec<String>)],
    parsed: &mut ParsedCifData,
) -> Result<()> {
    let mut columns = CifAtomSiteColumns::default();
    columns.group_idx = headers.iter().position(|h| h == "_atom_site.group_PDB");
    columns.id_idx = headers.iter().position(|h| h == "_atom_site.id");
    columns.name_idx = headers.iter().position(|h| h == "_atom_site.label_atom_id");
    if columns.name_idx.is_none() {
        columns.name_idx = headers.iter().position(|h| h == "_atom_site.auth_atom_id");
    }
    columns.resname_idx = headers.iter().position(|h| h == "_atom_site.label_comp_id");
    if columns.resname_idx.is_none() {
        columns.resname_idx = headers.iter().position(|h| h == "_atom_site.auth_comp_id");
    }
    columns.resnum_idx = headers.iter().position(|h| h == "_atom_site.label_seq_id");
    if columns.resnum_idx.is_none() {
        columns.resnum_idx = headers.iter().position(|h| h == "_atom_site.auth_seq_id");
    }
    columns.chain_idx = headers.iter().position(|h| h == "_atom_site.label_asym_id");
    if columns.chain_idx.is_none() {
        columns.chain_idx = headers.iter().position(|h| h == "_atom_site.auth_asym_id");
    }
    columns.x_idx = headers.iter().position(|h| h == "_atom_site.Cartn_x");
    columns.y_idx = headers.iter().position(|h| h == "_atom_site.Cartn_y");
    columns.z_idx = headers.iter().position(|h| h == "_atom_site.Cartn_z");

    let group_idx = columns.group_idx.context("Missing _atom_site.group_PDB in mmCIF loop")?;
    let asym_idx = columns.chain_idx.context("Missing _atom_site.label_asym_id in mmCIF loop")?;

    let mut serial = 1u32;
    for (ln, fields) in rows {
        let Some(group) = fields.get(group_idx).map(String::as_str) else { continue };
        if group != "ATOM" && group != "HETATM" {
            continue;
        }
        let Some(asym) = fields.get(asym_idx).cloned() else { continue };
        let atom = Atom::from_cif_atom_fields(
            &fields.iter().map(String::as_str).collect::<Vec<_>>(),
            &columns,
            serial,
        )
        .with_context(|| format!("Could not parse atom on line {ln}"))?;
        serial += 1;
        parsed.atoms_by_asym.entry(asym).or_default().push(atom);
    }
    Ok(())
}

fn parse_oper_list_loop_rows(
    headers: &[String],
    rows: &[(usize, Vec<String>)],
    parsed: &mut ParsedCifData,
) -> Result<()> {
    let id_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.id")
        .context("Missing _pdbx_struct_oper_list.id")?;
    let idx = |name: &str| headers.iter().position(|h| h == name).with_context(|| format!("Missing {name}"));
    let m11 = idx("_pdbx_struct_oper_list.matrix[1][1]")?;
    let m12 = idx("_pdbx_struct_oper_list.matrix[1][2]")?;
    let m13 = idx("_pdbx_struct_oper_list.matrix[1][3]")?;
    let v1 = idx("_pdbx_struct_oper_list.vector[1]")?;
    let m21 = idx("_pdbx_struct_oper_list.matrix[2][1]")?;
    let m22 = idx("_pdbx_struct_oper_list.matrix[2][2]")?;
    let m23 = idx("_pdbx_struct_oper_list.matrix[2][3]")?;
    let v2 = idx("_pdbx_struct_oper_list.vector[2]")?;
    let m31 = idx("_pdbx_struct_oper_list.matrix[3][1]")?;
    let m32 = idx("_pdbx_struct_oper_list.matrix[3][2]")?;
    let m33 = idx("_pdbx_struct_oper_list.matrix[3][3]")?;
    let v3 = idx("_pdbx_struct_oper_list.vector[3]")?;

    for (ln, fields) in rows {
        let Some(id) = fields.get(id_idx) else { continue };
        let parse = |idx: usize, field: &str| -> Result<f32> {
            fields
                .get(idx)
                .context(format!("Missing {field} on line {ln}"))?
                .parse()
                .with_context(|| format!("Could not parse {field} on line {ln}"))
        };
        let rot = [
            [parse(m11, "_pdbx_struct_oper_list.matrix[1][1]")?, parse(m12, "_pdbx_struct_oper_list.matrix[1][2]")?, parse(m13, "_pdbx_struct_oper_list.matrix[1][3]")?],
            [parse(m21, "_pdbx_struct_oper_list.matrix[2][1]")?, parse(m22, "_pdbx_struct_oper_list.matrix[2][2]")?, parse(m23, "_pdbx_struct_oper_list.matrix[2][3]")?],
            [parse(m31, "_pdbx_struct_oper_list.matrix[3][1]")?, parse(m32, "_pdbx_struct_oper_list.matrix[3][2]")?, parse(m33, "_pdbx_struct_oper_list.matrix[3][3]")?],
        ];
        let tr = Vec3::new(
            parse(v1, "_pdbx_struct_oper_list.vector[1]")?,
            parse(v2, "_pdbx_struct_oper_list.vector[2]")?,
            parse(v3, "_pdbx_struct_oper_list.vector[3]")?,
        ) / 10.0; // A -> nm
        parsed.operations.insert(id.clone(), Transform { rot, tr });
    }
    Ok(())
}

fn parse_assembly_gen_loop_rows(
    headers: &[String],
    rows: &[(usize, Vec<String>)],
    parsed: &mut ParsedCifData,
) -> Result<()> {
    let assembly_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_assembly_gen.assembly_id")
        .context("Missing _pdbx_struct_assembly_gen.assembly_id")?;
    let asym_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_assembly_gen.asym_id_list")
        .context("Missing _pdbx_struct_assembly_gen.asym_id_list")?;
    let oper_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_assembly_gen.oper_expression")
        .context("Missing _pdbx_struct_assembly_gen.oper_expression")?;

    for (_ln, fields) in rows {
        let (Some(assembly_id), Some(asym_list), Some(oper_expression)) =
            (fields.get(assembly_idx), fields.get(asym_idx), fields.get(oper_idx))
        else {
            continue;
        };
        let asym_ids = asym_list
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty() && s != "." && s != "?")
            .collect::<Vec<_>>();
        parsed.assembly_rules.push(AssemblyRule {
            assembly_id: assembly_id.clone(),
            asym_ids,
            oper_expression: oper_expression.clone(),
        });
    }
    Ok(())
}

fn build_cif_assembly(parsed: &ParsedCifData, assembly_id: &str) -> Result<Vec<Atom>> {
    let mut out = Vec::new();
    for rule in parsed.assembly_rules.iter().filter(|r| r.assembly_id == assembly_id) {
        let op_sequences = expand_oper_expression(&rule.oper_expression)?;
        for asym_id in &rule.asym_ids {
            let src = parsed
                .atoms_by_asym
                .get(asym_id)
                .with_context(|| format!("Referenced asym_id {asym_id} was not found"))?;
            for seq in &op_sequences {
                let tf = compose_operation_sequence(seq, &parsed.operations)?;
                for atom in src {
                    out.push(apply_transform(atom.clone(), tf));
                }
            }
        }
    }
    Ok(out)
}

fn compose_operation_sequence(seq: &[String], operations: &HashMap<String, Transform>) -> Result<Transform> {
    let mut composed = Transform::identity();
    for op_id in seq.iter().rev() {
        let op = operations.get(op_id).with_context(|| format!("Referenced operation id {op_id} was not found"))?;
        composed = compose(*op, composed);
    }
    Ok(composed)
}

fn compose(a: Transform, b: Transform) -> Transform {
    let mut rot = [[0.0f32; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            rot[r][c] =
                a.rot[r][0] * b.rot[0][c] + a.rot[r][1] * b.rot[1][c] + a.rot[r][2] * b.rot[2][c];
        }
    }
    let btr = apply_rot(a.rot, b.tr);
    Transform { rot, tr: btr + a.tr }
}

fn apply_transform(mut atom: Atom, tf: Transform) -> Atom {
    atom.pos = apply_rot(tf.rot, atom.pos) + tf.tr;
    atom
}

fn apply_rot(rot: [[f32; 3]; 3], v: Vec3) -> Vec3 {
    Vec3::new(
        rot[0][0] * v.x + rot[0][1] * v.y + rot[0][2] * v.z,
        rot[1][0] * v.x + rot[1][1] * v.y + rot[1][2] * v.z,
        rot[2][0] * v.x + rot[2][1] * v.y + rot[2][2] * v.z,
    )
}

fn expand_oper_expression(expr: &str) -> Result<Vec<Vec<String>>> {
    let cleaned: String = expr.chars().filter(|c| !c.is_whitespace()).collect();
    if cleaned.is_empty() {
        anyhow::bail!("Malformed oper_expression: empty");
    }
    if !cleaned.contains('(') {
        return Ok(vec![parse_oper_group(&cleaned)?]);
    }

    let mut groups: Vec<Vec<String>> = Vec::new();
    let bytes = cleaned.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] != b'(' {
            anyhow::bail!("Malformed oper_expression: {expr}");
        }
        i += 1;
        let start = i;
        while i < bytes.len() && bytes[i] != b')' {
            i += 1;
        }
        if i >= bytes.len() {
            anyhow::bail!("Malformed oper_expression: {expr}");
        }
        groups.push(parse_oper_group(&cleaned[start..i])?);
        i += 1;
    }

    let mut seqs: Vec<Vec<String>> = vec![Vec::new()];
    for group in groups {
        let mut next = Vec::new();
        for prefix in &seqs {
            for id in &group {
                let mut seq = prefix.clone();
                seq.push(id.clone());
                next.push(seq);
            }
        }
        seqs = next;
    }
    Ok(seqs)
}

fn parse_oper_group(group: &str) -> Result<Vec<String>> {
    if group.is_empty() {
        anyhow::bail!("Malformed oper_expression group");
    }
    let mut ids = Vec::new();
    for part in group.split(',') {
        if part.is_empty() {
            continue;
        }
        if let Some((start, end)) = part.split_once('-') {
            let s: i32 = start.parse().context("Malformed range start in oper_expression")?;
            let e: i32 = end.parse().context("Malformed range end in oper_expression")?;
            if s > e {
                anyhow::bail!("Malformed oper_expression range: {part}");
            }
            for v in s..=e {
                ids.push(v.to_string());
            }
        } else {
            ids.push(part.to_string());
        }
    }
    if ids.is_empty() {
        anyhow::bail!("Malformed oper_expression group");
    }
    Ok(ids)
}

#[derive(Default)]
struct CifAtomSiteColumns {
    group_idx: Option<usize>,
    id_idx: Option<usize>,
    name_idx: Option<usize>,
    resname_idx: Option<usize>,
    resnum_idx: Option<usize>,
    chain_idx: Option<usize>,
    x_idx: Option<usize>,
    y_idx: Option<usize>,
    z_idx: Option<usize>,
}
