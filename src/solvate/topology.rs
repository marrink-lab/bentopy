use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::{io::Read, path::PathBuf};

use anyhow::{Context as _, Result, bail};

pub fn determine_system_charge<P: AsRef<std::path::Path> + std::fmt::Debug>(
    topol_path: P,
) -> Result<i64> {
    let topol = {
        let mut file = std::fs::File::open(&topol_path)?;
        let mut buf = String::new();
        file.read_to_string(&mut buf)?;
        buf
    };

    // First, we want to find the [ system ] and [ molecules ] directives, because we need to take
    // inventory of the molecules we need to actually gather charge information for.
    type MoleculeName<'s> = &'s str;
    type MoleculeNumber = u64;
    let mut molecules: HashMap<MoleculeName, MoleculeNumber> = Default::default();
    // TODO: We're not actually implementing the ability to use backslashes to continue lines. That
    // would require a Cow.
    let mut lines = topol
        .lines()
        .map(strip_comments)
        .enumerate()
        .map(|(ln, line)| (ln + 1, line))
        .filter(|(_ln, line)| !line.is_empty());

    // "Using [ molecules ] without having used [ system ] before is meaningless and generates a
    // warning."
    skip_to_section(&mut lines, "system");
    // The next line is the system title.
    let (_ln, _title) = lines.next().context("Expected system title")?;

    // "After [ system ] the only allowed directive is [ molecules ]."
    skip_to_section(&mut lines, "molecules");
    // The next lines contain the list of molecules in the system.
    for (ln, line) in lines {
        // Parse a molecule.
        let Some((name, number)) = line.split_once(char::is_whitespace) else {
            bail!("Expected molecule declaration at line {ln}, found {line:?}");
        };
        let number: MoleculeNumber = number
            .trim_start()
            .parse()
            .with_context(|| format!("Expected an integer, found {number:?} at line {ln}"))?;
        *molecules.entry(name).or_default() += number;
    }

    // Now we go through the top file a second time, dilligently following all includes and
    // processing the relevant [ molecule_type ] sections.
    // That means that we are going to treat the root top file as if it were any normal itp file,
    // and we follow every include in a depth-first manner.
    let moleculetypes = Context::parse_molecule_types(topol_path.as_ref())
        .with_context(|| format!("Problem while parsing main topology file {topol_path:?}"))?;

    // Finally, the accounting. We now know how many of each moleculetype we have, and what the
    // charge of each individual instance of that moleculetype is. Now we compute the total charge.
    let mut total_charge = 0.0;
    for (name, number) in molecules {
        let Some(charge): Option<f64> =
            moleculetypes.iter().find_map(|m| if m.name == name { Some(m.charge) } else { None })
        else {
            bail!("No moleculetype was defined for {name:?}");
        };

        // Report any fractional charges of individual structures as warnings. This makes it much
        // easier to identify issues since we have a hard error on fractional charges below.
        if is_sufficiently_integer_charge(charge) {
            eprintln!("WARNING: Molecule {name:?} has a fractional charge of {charge}.");
        }

        total_charge += number as f64 * charge;
    }

    if is_sufficiently_integer_charge(total_charge) {
        bail!("The system charge should be an integer value, but it is {total_charge}");
    }
    let total_charge = total_charge.round() as i64;

    Ok(total_charge)
}

/// Trim any surrounding whitespace and strip comments.
fn strip_comments(line: &str) -> &str {
    let line = line.trim_start();
    match line.split_once(';') {
        Some((content, _comment)) => content,
        None => line,
    }
    .trim_end()
}

fn parse_directive(line: &str) -> Option<&str> {
    line.strip_prefix('[').and_then(|rest| rest.strip_suffix(']')).map(str::trim)
}

fn parse_path_str<'s>(s: &'s str, ln: usize, line: &'s str) -> Result<&'s str> {
    let front = s.trim_start().strip_prefix('"').with_context(|| {
        format!("Expected the start of a quoted file path at line {ln}, found {line:?}")
    })?;
    front.strip_suffix('"').with_context(|| {
        format!("Expected the end of a quoted file path at line {ln}, found {line:?}")
    })
}

fn skip_while<'s>(
    lines: &mut impl Iterator<Item = (usize, &'s str)>,
    condition: impl Fn(&'s str) -> bool,
) {
    while let Some((_ln, line)) = lines.next()
        && condition(line)
    {}
}

fn skip_to_section<'s>(
    lines: &mut impl Iterator<Item = (usize, &'s str)>,
    directive: &'static str,
) {
    skip_while(lines, |line| parse_directive(line) != Some(directive))
}

fn is_sufficiently_integer_charge(charge: f64) -> bool {
    /// The maximum allowable deviation from an integer for the total charge of a system.
    ///
    /// This value was picked in an informed but ultimately arbitrary way. Open to more deeply reasoned
    /// alternatives.
    const FRACTIONAL_CHARGE_THRESHOLD: f64 = 0.01;
    f64::fract(charge.round() - charge).abs() > FRACTIONAL_CHARGE_THRESHOLD
}

struct Molecule {
    name: String,
    charge: f64,
}

impl Molecule {
    fn new(name: String) -> Self {
        Self { name, charge: Default::default() }
    }
}

/// Parser state.
#[derive(Clone, Copy, Default)]
enum State {
    /// Not in a (relevant) section.
    #[default]
    None,
    /// In an `[ atomtypes ]` section.
    AtomTypes,
    /// In a `[ moleculetype ]` section.
    MoleculeType,
    /// In an `[ atoms ]` section.
    Atoms,
}

#[derive(Default)]
struct Context {
    // TODO: Small string?
    /// The set of defined parameters.
    defined: HashSet<String>,
    /// Handle nested ifdefs, tracking if they are active (`true`) or inactive (`false`) for
    /// managing correct handling of `#else` directives.
    stack: Vec<bool>,
    // TODO: Small string?
    atoms: HashMap<String, f64>,
    molecules: Vec<Molecule>,
    current: Option<Molecule>,
    state: State,
}

impl Context {
    fn is_active(&self) -> bool {
        self.stack.last().copied().unwrap_or(true)
    }

    fn parse(&mut self, path: &Path) -> Result<()> {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        let parent = path.parent().unwrap_or(std::path::Path::new("."));

        for (ln, line) in reader.lines().enumerate() {
            let ln = ln + 1;
            let line = line?;
            let line = strip_comments(&line);
            if line.is_empty() {
                continue;
            }

            // Move this to a method and slap a context around it for line numbers etc.
            if let Some(hash_directive_line) = line.strip_prefix('#') {
                let (directive, rest) = hash_directive_line
                    .split_once(char::is_whitespace)
                    .unwrap_or((hash_directive_line, ""));
                match directive {
                    // Immediately follow any #include directives.
                    "include" => {
                        let relative_path = parse_path_str(rest, ln, line)?;
                        let mut full_path = PathBuf::from(parent);
                        full_path.push(relative_path);
                        self.parse(&full_path)
                            .with_context(|| format!("Could not parse {full_path:?}"))?;
                    }

                    // Change the conditionals stack.
                    "ifdef" => {
                        let symbol =
                            rest.split_whitespace().next().context("Expected a definition name")?;
                        // We go a level deeper.
                        let current = self.defined.contains(symbol);
                        let parent_active = self.is_active();
                        self.stack.push(current && parent_active);
                    }
                    "ifndef" => {
                        let symbol =
                            rest.split_whitespace().next().context("Expected a definition name")?;
                        // We go a level deeper.
                        let current = !self.defined.contains(symbol);
                        let parent_active = self.is_active();
                        self.stack.push(current && parent_active);
                    }
                    "else" => {
                        let prev = self.stack.pop().context(
                            "Encountered #else directive without preceding #ifdef or #ifndef",
                        )?;
                        // Invert our activity state in the else branch.
                        let parent_active = self.is_active();
                        self.stack.push(parent_active && !prev);
                    }
                    "endif" => {
                        self.stack.pop().context("Popped too far!")?; // Pop down a level.
                    }

                    // Important to not consider any #define or #undef directives when we are
                    // not active.
                    "define" if self.is_active() => {
                        let symbol =
                            rest.split_whitespace().next().context("Expected a definition name")?;
                        self.defined.insert(symbol.to_string());
                    }
                    "undef" if self.is_active() => {
                        let symbol =
                            rest.split_whitespace().next().context("Expected a definition name")?;
                        self.defined.insert(symbol.to_string());
                    }

                    unknown => {
                        eprintln!("UNKNOWN DIRECTIVE: #{unknown} at {ln} in {path:?}:\n{line:?}");
                    }
                }

                continue;
            }

            // Don't read any lines if we're not active.
            if !self.is_active() {
                continue;
            }

            // Consider if we are in a new section.
            if let Some(directive) = parse_directive(line) {
                match directive.to_lowercase().as_str() {
                    "atomtypes" => self.state = State::AtomTypes,
                    "moleculetype" => {
                        self.wrap_up();
                        self.state = State::MoleculeType;
                    }
                    "atoms" => self.state = State::Atoms,
                    _ => self.state = State::None,
                }
                continue;
            }

            // Consume section lines.
            match self.state {
                State::None => {}
                State::AtomTypes => self.parse_atom_parameters(line).with_context(|| {
                    format!("Could not parse atom parameters line at {ln} in {path:?}: {line:?}")
                })?,
                State::MoleculeType => {
                    // TODO: We're leaving a small performance improvement on the table here. We
                    // can check whether the name is actually relevant for our final tally, and
                    // depending on that choose to read or not read the molecule definition.
                    if let Some(name) = line.split_whitespace().next() {
                        self.current = Some(Molecule::new(name.to_string()));
                    }
                }
                State::Atoms => self.parse_molecule_atom(line).with_context(|| {
                    let name = self.current.as_ref().map(|curr| curr.name.as_str()).unwrap_or("?");
                    format!(
                        "Could not parse atom line at {ln} in {path:?} \
                                    for molecule {name:?}: {line:?}"
                    )
                })?,
            }
        }

        Ok(())
    }

    fn parse_atom_parameters(&mut self, line: &str) -> Result<()> {
        // Because of cursed GROMACS shit, we need to test around to see where the charge field
        // is. We need to check for the field being a 'Particle Type'.
        // https://gitlab.com/gromacs/gromacs/-/blob/6f48ea2adfb7263c99f443f8d54b5b0f1da33aa6/src/gromacs/gmxpreprocess/toppush.cpp#L351-379
        /* Comments on optional fields in the atomtypes section:
         *
         * The force field format is getting a bit old. For OPLS-AA we needed
         * to add a special bonded atomtype, and for Gerrit Groenhofs QM/MM stuff
         * we also needed the atomic numbers.
         * To avoid making all old or user-generated force fields unusable we
         * have introduced both these quantities as optional columns, and do some
         * acrobatics to check whether they are present or not.
         * This will all look much nicer when we switch to XML... sigh.
         *
         * Field 0 (mandatory) is the nonbonded type name. (string)
         * Field 1 (optional)  is the bonded type (string)
         * Field 2 (optional)  is the atomic number (int)
         * Field 3 (mandatory) is the mass (numerical)
         * Field 4 (mandatory) is the charge (numerical)
         * Field 5 (mandatory) is the particle type (single character)
         * This is followed by a number of nonbonded parameters.
         *
         * The safest way to identify the format is the particle type field.
         *
         * So, here is what we do:
         *
         * A. Read in the first six fields as strings
         * B. If field 3 (starting from 0) is a single char, we have neither
         *    bonded_type or atomic numbers.
         * C. If field 5 is a single char we have both.
         * D. If field 4 is a single char we check field 1. If this begins with
         *    an alphabetical character we have bonded types, otherwise atomic numbers.
         */

        // A. Read in the first six fields as strings
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 6 {
            bail!("Too few fields in atom type definition");
        }

        let field_is_single_char = |idx: usize| {
            fields[idx].len() == 1 && fields[idx].chars().next().unwrap().is_ascii_alphabetic()
        };
        let charge_idx = if field_is_single_char(5) {
            // C. If field 5 is a single char we have both.
            4
        } else if field_is_single_char(3) {
            // B. If field 3 (starting from 0) is a single char, we have neither
            //    bonded_type or atomic numbers.
            2
        } else if field_is_single_char(4) {
            // D. If field 4 is a single char we check field 1. If this begins with
            //    an alphabetical character we have bonded types, otherwise atomic numbers.
            3
        } else {
            bail!("The particle type for this atom line cannot be reliably determined")
        };

        // TODO: Small string optimization?
        let atom_type_name = fields[0].to_string();
        let charge = fields[charge_idx];
        let charge = charge.parse().with_context(|| {
            format!("could not parse charge value {charge:?} for {atom_type_name} as a float")
        })?;

        self.atoms.insert(atom_type_name.to_string(), charge);
        Ok(())
    }

    fn parse_molecule_atom(&mut self, line: &str) -> Result<()> {
        let mut fields = line.split_whitespace();
        let Some(current) = self.current.as_mut() else {
            bail!("Cannot read an atoms section without a prior [ moleculetype ] directive")
        };

        // Fields: nr, atom type, residue number, residue name, atom name, charge group number, q(e), m(u).
        let atom_type_name = fields.nth(1).context("Expected atom type name")?;
        // If the moleculetype atom contains a charge, that overrides the charge set in the
        // forcefield parameters for that atomtype.
        let charge = if let Some(charge) = fields.nth(4) {
            charge
                .parse()
                .with_context(|| format!("Could not parse charge value {charge:?} in atom line"))?
        } else if let Some(&charge) = self.atoms.get(atom_type_name) {
            charge
        } else {
            // I think this would be considered an error in the eyes of GROMACS as well.
            bail!(
                "Atom type {atom_type_name} was not previously defined \
                            and has no charge value in this atom line"
            )
        };
        current.charge += charge;
        Ok(())
    }

    fn wrap_up(&mut self) {
        if let Some(current) = self.current.take() {
            self.molecules.push(current);
        }
    }

    fn parse_molecule_types(topol_path: &Path) -> Result<Vec<Molecule>> {
        let mut context = Context::default();
        context.parse(topol_path)?;
        context.wrap_up();
        Ok(context.molecules)
    }
}
