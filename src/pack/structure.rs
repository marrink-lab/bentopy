use std::io::{BufRead, BufReader};
use std::num::NonZeroUsize;
use std::collections::HashMap;

use eightyseven::reader::{ParseList, ReadGro};
use eightyseven::structure::{AtomName, AtomNum, ResName, ResNum};
use eightyseven::writer::{WriteGro, format_atom_line};
use glam::Vec3;

pub type Atom = Vec3;

/// Load a [`Structure`] from a structure file.
///
/// This function will return an error if the structure is empty. Downstream functions may assume
/// that a [`Structure`] has at least one atom.
pub fn load_molecule<P: AsRef<std::path::Path> + std::fmt::Debug>(
    path: P,
) -> anyhow::Result<Structure> {
    use anyhow::{Context, bail};

    let file = std::fs::File::open(&path)?;

    let structure = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gro") => Structure::read_from_file(file)
            .with_context(|| format!("Failed to parse gro file {path:?}"))?,
        Some("pdb") => Structure::read_from_pdb_file(file)
            .with_context(|| format!("Failed to parse PDB file {path:?}"))?,
        Some("cif") | Some("mmcif") => Structure::read_from_cif_file(file)
            .with_context(|| format!("Failed to parse mmCIF file {path:?}"))?,
        None | Some(_) => {
            eprintln!("WARNING: Assuming {path:?} is a PDB file.");
            Structure::read_from_pdb_file(file)
                .with_context(|| format!("Failed to parse the file {path:?} as PDB"))?
        }
    };

    if structure.natoms() == 0 {
        bail!("Structure from {path:?} contains no atoms")
    }

    Ok(structure)
}

/// A structure type that stores its atoms as simple positions.
///
/// Invariant: A `Structure` is always centered, such that its geometric center lies at the origin.
///
/// Invariant: A `Structure` has at least one atom.
pub struct Structure {
    atoms: Vec<Atom>,
}

impl ReadGro<Atom> for Structure {
    const PARSE_LIST: ParseList = ParseList {
        resnum: false,
        resname: false,
        atomname: false,
        atomnum: false,
        position: true,
        velocity: false,
    };

    fn build_atom(
        _resnum: Option<ResNum>,
        _resname: Option<ResName>,
        _atomname: Option<AtomName>,
        _atomnum: Option<AtomNum>,
        position: Option<[f32; 3]>,
        _velocity: Option<[f32; 3]>,
    ) -> Atom {
        // We can safely unwrap because this is the only value we expect.
        Atom::from_array(position.unwrap())
    }

    fn build_structure(
        _title: String,
        atoms: Vec<Atom>,
        _boxvecs: eightyseven::structure::BoxVecs,
    ) -> Self {
        let mut structure = Self { atoms };
        structure.translate_to_center();
        structure
    }
}

#[non_exhaustive]
#[derive(Debug)]
pub enum ParsePdbError {
    IOError(std::io::Error),
    BadPositionValue {
        err: std::num::ParseFloatError,
        /// Line number.
        ln: NonZeroUsize,
    },
}

impl std::error::Error for ParsePdbError {}

impl std::fmt::Display for ParsePdbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IOError(err) => write!(f, "IO error: {err}"),
            Self::BadPositionValue { err, ln } => {
                write!(f, "Could not parse record on line {ln}: {err}")
            }
        }
    }
}

impl From<std::io::Error> for ParsePdbError {
    fn from(err: std::io::Error) -> Self {
        Self::IOError(err)
    }
}

#[non_exhaustive]
#[derive(Debug)]
pub enum ParseCifError {
    IOError(std::io::Error),
    MissingAtomSiteField(&'static str),
    MissingOperListField(&'static str),
    MissingAssemblyGenField(&'static str),
    MissingOperListValue {
        field: &'static str,
        ln: NonZeroUsize,
    },
    MissingOperationId(String),
    MissingAsymId(String),
    MalformedOperExpression(String),
    BadPositionValue {
        err: std::num::ParseFloatError,
        /// Line number.
        ln: NonZeroUsize,
    },
}

impl std::error::Error for ParseCifError {}

impl std::fmt::Display for ParseCifError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IOError(err) => write!(f, "IO error: {err}"),
            Self::MissingAtomSiteField(field) => {
                write!(f, "Could not find required _atom_site field {field}")
            }
            Self::MissingOperListField(field) => {
                write!(f, "Could not find required _pdbx_struct_oper_list field {field}")
            }
            Self::MissingAssemblyGenField(field) => {
                write!(f, "Could not find required _pdbx_struct_assembly_gen field {field}")
            }
            Self::MissingOperListValue { field, ln } => {
                write!(f, "Could not parse _pdbx_struct_oper_list field {field} on line {ln}")
            }
            Self::MissingOperationId(id) => {
                write!(f, "Could not find referenced operation id {id}")
            }
            Self::MissingAsymId(id) => {
                write!(f, "Could not find referenced asym_id {id}")
            }
            Self::MalformedOperExpression(expr) => {
                write!(f, "Malformed oper_expression {expr}")
            }
            Self::BadPositionValue { err, ln } => {
                write!(f, "Could not parse _atom_site record on line {ln}: {err}")
            }
        }
    }
}

impl From<std::io::Error> for ParseCifError {
    fn from(err: std::io::Error) -> Self {
        Self::IOError(err)
    }
}

impl Structure {
    pub fn read_from_pdb_file(file: std::fs::File) -> Result<Structure, ParsePdbError> {
        let reader = BufReader::new(file);

        let mut atoms = Vec::new();
        for (ln, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let ln = NonZeroUsize::new(ln + 1).unwrap(); // We know ln + 1 >= 1 because ln >= 0.
                let bad_position = |err| ParsePdbError::BadPositionValue { err, ln };
                let x = line[30..38].trim().parse().map_err(bad_position)?;
                let y = line[38..46].trim().parse().map_err(bad_position)?;
                let z = line[46..54].trim().parse().map_err(bad_position)?;
                let atom = Atom::new(x, y, z) / 10.0; // Convert from Å to nm.
                atoms.push(atom);
            }
        }

        let mut structure = Self { atoms };
        structure.translate_to_center();
        Ok(structure)
    }

    pub fn read_from_cif_file(file: std::fs::File) -> Result<Structure, ParseCifError> {
        let parsed = parse_cif(file)?;
        // Prefer biological assembly "1" when assembly metadata is present.
        if parsed.operations.is_empty() || parsed.assembly_gen.is_empty() {
            build_asymmetric_unit(parsed)
        } else {
            build_with_assembly(parsed, "1")
        }
    }

    /// Translate this [`Structure`] such that it's geometric center lies at the origin.
    fn translate_to_center(&mut self) {
        // Invariant: A Structure has at least one atom.
        let center = self.atoms().sum::<Atom>() / self.natoms() as f32;
        for pos in &mut self.atoms {
            *pos -= center;
        }
    }

    /// Calculate the moment of inertia for this [`Structure`].
    ///
    /// Invariant: Assumes that the structure is centered.
    ///
    ///     I = Σ(m_i * r_i²)
    pub fn moment_of_inertia(&self) -> f32 {
        self.atoms().map(|atom| atom.length_squared()).sum()
    }

    /// Returns the radius of the point that is farthest from the structure geometric center.
    ///
    /// Invariant: Assumes that the structure is centered.
    pub fn bounding_sphere(&self) -> f32 {
        self.atoms().map(|atom| atom.length()).max_by(f32::total_cmp).unwrap() // Invariant: A structure has at least one atom.
    }
}

impl<'atoms> WriteGro<'atoms, Atom> for Structure {
    fn title(&self) -> String {
        "debug".to_string()
    }

    fn natoms(&self) -> usize {
        self.atoms.len()
    }

    fn atoms(&'atoms self) -> impl Iterator<Item = &'atoms Atom> {
        self.atoms.iter()
    }

    fn boxvecs(&self) -> String {
        "400.0 400.0 400.0".to_string()
    }

    fn format_atom_line(atom: &Atom) -> String {
        format_atom_line(1, "DUMMY", "DUMMY", 2, atom.to_array(), None)
    }
}

fn build_asymmetric_unit(parsed: ParsedCif) -> Result<Structure, ParseCifError> {
    let mut atoms = Vec::new();
    for chain_atoms in parsed.atoms_by_asym.values() {
        atoms.extend_from_slice(chain_atoms);
    }

    if atoms.is_empty() {
        return Ok(Structure { atoms });
    }

    let mut structure = Structure { atoms };
    structure.translate_to_center();
    Ok(structure)
}

fn build_with_assembly(parsed: ParsedCif, assembly_id: &str) -> Result<Structure, ParseCifError> {
    let mut atoms = Vec::new();
    for rule in parsed.assembly_gen.iter().filter(|r| r.assembly_id == assembly_id) {
        let op_sequences = expand_oper_expression(&rule.oper_expression)?;
        for asym_id in &rule.asym_ids {
            let source_atoms = parsed
                .atoms_by_asym
                .get(asym_id)
                .ok_or_else(|| ParseCifError::MissingAsymId(asym_id.clone()))?;
            for seq in &op_sequences {
                let transform = compose_operation_sequence(seq, &parsed.operations)?;
                for atom in source_atoms {
                    atoms.push(apply_transform(*atom, transform));
                }
            }
        }
    }

    if atoms.is_empty() {
        return Ok(Structure { atoms });
    }

    let mut structure = Structure { atoms };
    structure.translate_to_center();
    Ok(structure)
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

struct AssemblyGenRule {
    assembly_id: String,
    asym_ids: Vec<String>,
    oper_expression: String,
}

struct ParsedCif {
    atoms_by_asym: HashMap<String, Vec<Atom>>,
    operations: HashMap<String, Transform>,
    assembly_gen: Vec<AssemblyGenRule>,
}

fn parse_cif(file: std::fs::File) -> Result<ParsedCif, ParseCifError> {
    let reader = BufReader::new(file);
    let lines: Vec<(NonZeroUsize, String)> = reader
        .lines()
        .enumerate()
        .map(|(i, line)| {
            let ln = NonZeroUsize::new(i + 1).unwrap();
            line.map(|l| (ln, l))
        })
        .collect::<Result<_, _>>()?;

    let mut atoms_by_asym: HashMap<String, Vec<Atom>> = HashMap::new();
    let mut operations = HashMap::new();
    let mut assembly_gen = Vec::new();

    let mut i = 0usize;
    while i < lines.len() {
        let trimmed = lines[i].1.trim();
        if trimmed != "loop_" {
            i += 1;
            continue;
        }

        let (headers, rows, next_i) = parse_loop_rows(&lines, i + 1);
        i = next_i;
        if headers.is_empty() {
            continue;
        }

        let first = headers[0].as_str();
        if first.starts_with("_atom_site.") {
            parse_atom_site_loop(&headers, &rows, &mut atoms_by_asym)?;
        } else if first.starts_with("_pdbx_struct_oper_list.") {
            parse_oper_list_loop(&headers, &rows, &mut operations)?;
        } else if first.starts_with("_pdbx_struct_assembly_gen.") {
            parse_assembly_gen_loop(&headers, &rows, &mut assembly_gen)?;
        }
    }

    Ok(ParsedCif { atoms_by_asym, operations, assembly_gen })
}

fn parse_loop_rows(
    lines: &[(NonZeroUsize, String)],
    mut i: usize,
) -> (Vec<String>, Vec<(NonZeroUsize, Vec<String>)>, usize) {
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
    let mut row_start_ln: Option<NonZeroUsize> = None;
    while i < lines.len() {
        let (ln, raw) = &lines[i];
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
                row_start_ln = Some(*ln);
            }
            buffer.extend(tokens);
            while ncols > 0 && buffer.len() >= ncols {
                let fields: Vec<String> = buffer.drain(0..ncols).collect();
                rows.push((row_start_ln.unwrap_or(*ln), fields));
                row_start_ln = if buffer.is_empty() { None } else { Some(*ln) };
            }
        }
        i += 1;
    }

    (headers, rows, i)
}

fn parse_atom_site_loop(
    headers: &[String],
    rows: &[(NonZeroUsize, Vec<String>)],
    atoms_by_asym: &mut HashMap<String, Vec<Atom>>,
) -> Result<(), ParseCifError> {
    let group_idx = headers
        .iter()
        .position(|h| h == "_atom_site.group_PDB")
        .ok_or(ParseCifError::MissingAtomSiteField("_atom_site.group_PDB"))?;
    let asym_idx = headers
        .iter()
        .position(|h| h == "_atom_site.label_asym_id")
        .ok_or(ParseCifError::MissingAtomSiteField("_atom_site.label_asym_id"))?;
    let x_idx = headers
        .iter()
        .position(|h| h == "_atom_site.Cartn_x")
        .ok_or(ParseCifError::MissingAtomSiteField("_atom_site.Cartn_x"))?;
    let y_idx = headers
        .iter()
        .position(|h| h == "_atom_site.Cartn_y")
        .ok_or(ParseCifError::MissingAtomSiteField("_atom_site.Cartn_y"))?;
    let z_idx = headers
        .iter()
        .position(|h| h == "_atom_site.Cartn_z")
        .ok_or(ParseCifError::MissingAtomSiteField("_atom_site.Cartn_z"))?;

    for (ln, fields) in rows {
        let Some(group) = fields.get(group_idx).map(String::as_str) else { continue };
        if group != "ATOM" && group != "HETATM" {
            continue;
        }
        let Some(asym) = fields.get(asym_idx) else { continue };
        if asym == "." || asym == "?" {
            continue;
        }
        let bad = |err| ParseCifError::BadPositionValue { err, ln: *ln };
        let Some(xs) = fields.get(x_idx) else { continue };
        let Some(ys) = fields.get(y_idx) else { continue };
        let Some(zs) = fields.get(z_idx) else { continue };
        let x: f32 = xs.parse().map_err(bad)?;
        let y: f32 = ys.parse().map_err(bad)?;
        let z: f32 = zs.parse().map_err(bad)?;
        atoms_by_asym
            .entry(asym.clone())
            .or_default()
            .push(Vec3::new(x, y, z) / 10.0); // A -> nm
    }

    Ok(())
}

fn parse_oper_list_loop(
    headers: &[String],
    rows: &[(NonZeroUsize, Vec<String>)],
    operations: &mut HashMap<String, Transform>,
) -> Result<(), ParseCifError> {
    let id_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.id")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.id"))?;
    let m11 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[1][1]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[1][1]"))?;
    let m12 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[1][2]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[1][2]"))?;
    let m13 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[1][3]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[1][3]"))?;
    let v1 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.vector[1]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.vector[1]"))?;
    let m21 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[2][1]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[2][1]"))?;
    let m22 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[2][2]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[2][2]"))?;
    let m23 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[2][3]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[2][3]"))?;
    let v2 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.vector[2]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.vector[2]"))?;
    let m31 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[3][1]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[3][1]"))?;
    let m32 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[3][2]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[3][2]"))?;
    let m33 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.matrix[3][3]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.matrix[3][3]"))?;
    let v3 = headers
        .iter()
        .position(|h| h == "_pdbx_struct_oper_list.vector[3]")
        .ok_or(ParseCifError::MissingOperListField("_pdbx_struct_oper_list.vector[3]"))?;

    for (ln, fields) in rows {
        let Some(id) = fields.get(id_idx) else { continue };
        let parse =
            |idx: usize, field: &'static str| -> Result<f32, ParseCifError> {
                match fields.get(idx) {
                    Some(v) => v.parse().map_err(|_| ParseCifError::MissingOperListValue {
                        field,
                        ln: *ln,
                    }),
                    None => Err(ParseCifError::MissingOperListValue { field, ln: *ln }),
                }
            };
        let rot = [
            [
                parse(m11, "_pdbx_struct_oper_list.matrix[1][1]")?,
                parse(m12, "_pdbx_struct_oper_list.matrix[1][2]")?,
                parse(m13, "_pdbx_struct_oper_list.matrix[1][3]")?,
            ],
            [
                parse(m21, "_pdbx_struct_oper_list.matrix[2][1]")?,
                parse(m22, "_pdbx_struct_oper_list.matrix[2][2]")?,
                parse(m23, "_pdbx_struct_oper_list.matrix[2][3]")?,
            ],
            [
                parse(m31, "_pdbx_struct_oper_list.matrix[3][1]")?,
                parse(m32, "_pdbx_struct_oper_list.matrix[3][2]")?,
                parse(m33, "_pdbx_struct_oper_list.matrix[3][3]")?,
            ],
        ];
        let tr = Vec3::new(
            parse(v1, "_pdbx_struct_oper_list.vector[1]")?,
            parse(v2, "_pdbx_struct_oper_list.vector[2]")?,
            parse(v3, "_pdbx_struct_oper_list.vector[3]")?,
        ) / 10.0; // A -> nm
        operations.insert(id.clone(), Transform { rot, tr });
    }

    Ok(())
}

fn parse_assembly_gen_loop(
    headers: &[String],
    rows: &[(NonZeroUsize, Vec<String>)],
    assembly_gen: &mut Vec<AssemblyGenRule>,
) -> Result<(), ParseCifError> {
    let assembly_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_assembly_gen.assembly_id")
        .ok_or(ParseCifError::MissingAssemblyGenField(
            "_pdbx_struct_assembly_gen.assembly_id",
        ))?;
    let asym_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_assembly_gen.asym_id_list")
        .ok_or(ParseCifError::MissingAssemblyGenField(
            "_pdbx_struct_assembly_gen.asym_id_list",
        ))?;
    let oper_idx = headers
        .iter()
        .position(|h| h == "_pdbx_struct_assembly_gen.oper_expression")
        .ok_or(ParseCifError::MissingAssemblyGenField(
            "_pdbx_struct_assembly_gen.oper_expression",
        ))?;

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

        assembly_gen.push(AssemblyGenRule {
            assembly_id: assembly_id.clone(),
            asym_ids,
            oper_expression: oper_expression.clone(),
        });
    }

    Ok(())
}

fn compose_operation_sequence(
    seq: &[String],
    operations: &HashMap<String, Transform>,
) -> Result<Transform, ParseCifError> {
    let mut composed = Transform::identity();
    for op_id in seq.iter().rev() {
        let op = operations
            .get(op_id)
            .ok_or_else(|| ParseCifError::MissingOperationId(op_id.clone()))?;
        composed = compose(*op, composed);
    }
    Ok(composed)
}

fn compose(a: Transform, b: Transform) -> Transform {
    // Compose as a ∘ b (apply b, then a)
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

fn apply_transform(atom: Atom, tf: Transform) -> Atom {
    apply_rot(tf.rot, atom) + tf.tr
}

fn apply_rot(rot: [[f32; 3]; 3], v: Vec3) -> Vec3 {
    Vec3::new(
        rot[0][0] * v.x + rot[0][1] * v.y + rot[0][2] * v.z,
        rot[1][0] * v.x + rot[1][1] * v.y + rot[1][2] * v.z,
        rot[2][0] * v.x + rot[2][1] * v.y + rot[2][2] * v.z,
    )
}

fn expand_oper_expression(expr: &str) -> Result<Vec<Vec<String>>, ParseCifError> {
    let cleaned: String = expr.chars().filter(|c| !c.is_whitespace()).collect();
    if cleaned.is_empty() {
        return Err(ParseCifError::MalformedOperExpression(expr.to_string()));
    }

    if !cleaned.contains('(') {
        return Ok(vec![parse_oper_group(&cleaned, expr)?]);
    }

    let mut groups: Vec<Vec<String>> = Vec::new();
    let bytes = cleaned.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] != b'(' {
            return Err(ParseCifError::MalformedOperExpression(expr.to_string()));
        }
        i += 1;
        let start = i;
        while i < bytes.len() && bytes[i] != b')' {
            i += 1;
        }
        if i >= bytes.len() {
            return Err(ParseCifError::MalformedOperExpression(expr.to_string()));
        }
        let inner = &cleaned[start..i];
        groups.push(parse_oper_group(inner, expr)?);
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

fn parse_oper_group(group: &str, original_expr: &str) -> Result<Vec<String>, ParseCifError> {
    if group.is_empty() {
        return Err(ParseCifError::MalformedOperExpression(original_expr.to_string()));
    }

    let mut ids = Vec::new();
    for part in group.split(',') {
        if part.is_empty() {
            continue;
        }
        if let Some((start, end)) = part.split_once('-') {
            let s: i32 = start
                .parse()
                .map_err(|_| ParseCifError::MalformedOperExpression(original_expr.to_string()))?;
            let e: i32 = end
                .parse()
                .map_err(|_| ParseCifError::MalformedOperExpression(original_expr.to_string()))?;
            if s > e {
                return Err(ParseCifError::MalformedOperExpression(original_expr.to_string()));
            }
            for v in s..=e {
                ids.push(v.to_string());
            }
        } else {
            ids.push(part.to_string());
        }
    }
    if ids.is_empty() {
        return Err(ParseCifError::MalformedOperExpression(original_expr.to_string()));
    }
    Ok(ids)
}
