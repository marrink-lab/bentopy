use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;
use std::{io::Read, path::PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};
use eightyseven::structure::ResName;

/// Solvate.
#[derive(Debug, Parser)]
#[command(about, version = bentopy::core::version::VERSION)]
pub struct Args {
    /// Structure input path.
    #[arg(short, long)]
    pub input: PathBuf,
    /// Output path.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Lowest allowable distance between solvent and structure beads (nm).
    ///
    /// This is the minimum allowable center-to-center distance when checking for a collision
    /// between a solvent bead and a structure bead.
    ///
    /// For the solvent-solvent cutoff distance, see `--solvent-cutoff`.
    ///
    /// For Martini solvation, the default cutoff is 0.43 nm. For atomistic solvation, the default
    /// cutoff is 0.28.
    #[arg(long)]
    pub cutoff: Option<f32>,

    /// Lowest allowable distance between solvent beads (nm).
    ///
    /// This is the minimum allowable center-to-center distance between solvent beads when cutting
    /// the sides to fit in the specified output structure box.
    ///
    /// For Martini solvation, the default cutoff is 0.21 nm. For atomistic solvation, the default
    /// cutoff is 0.21.
    #[arg(long)]
    pub solvent_cutoff: Option<f32>,

    /// List of resnames to ignore when checking against structure-solvent collisions.
    #[arg(long, value_delimiter = ',')]
    pub ignore: Vec<ResName>,

    /// The type of water written to the output file.
    #[arg(long, value_enum, default_value_t)]
    pub water_type: WaterType,

    /// Center the structure in the new box.
    ///
    /// Note that this option only has an effect if the boundary mode is set to `grow`.
    #[arg(short, long)]
    pub center: bool,

    /// Set how the boundaries of the box will be treated.
    #[arg(short, long, value_enum, default_value_t)]
    pub boundary_mode: BoundaryMode,

    /// Set how periodic images of the input structure are considered.
    ///
    /// Note that only orthorhombic (cuboid) periodicity is considered, currently.
    // TODO: Consider whether supporting other forms of periodicity is worthwhile.
    #[arg(short, long, value_enum, default_value_t)]
    pub periodic_mode: PeriodicMode,

    /// Define substitutes for solvent residues, such as ions.
    #[arg(short, long = "substitute")]
    pub substitutes: Vec<SubstituteConfig>,

    /// Set the charge to neutralize with additional ions.
    ///
    /// In essence, this is a shorthand for explicitly providing ions to compensate the charge as
    /// substitutes. This can be helpful for automation purposes.
    ///
    /// When you pass a charge value of 'auto', the residual charge will be determined on a
    /// best-effort basis by reading the topology file, if specified.
    ///
    /// By default, NA (positive) and CL (negative) ions are used, but a different
    /// positive-negative substitution pair can be specified by prepending a colon followed by the
    /// positive and negative substitute name separated by one comma.
    ///
    ///     <charge>:<positive ion name>,<negative ion name>
    ///
    /// So, by default, some `<charge>` is interpreted as
    ///
    ///     <charge>:NA,CL
    #[arg(long, allow_hyphen_values = true)]
    pub charge: Option<ChargeConfig>,

    /// Combine substitutes with identical names into one block.
    #[arg(long, default_value_t)]
    pub no_combine_substitutes: bool,

    /// Set whether and in what way substitutes are sorted.
    #[arg(long, value_enum, default_value_t)]
    pub sort_substitutes: SortBehavior,

    /// Random number generator seed for solvent substitution, such as ion placement.
    #[arg(long)]
    pub seed: Option<u64>,

    /// If the solvent template contains velocity information, write these velocities to the output
    /// file.
    #[arg(long, default_value_t)]
    pub write_velocities: bool,

    /// Append solvation topology lines to a path.
    #[arg(short = 't', long)]
    pub append_topol: Option<PathBuf>,

    #[arg(long)]
    pub no_write_parallel: bool,
    /// The suggested number of atoms to format at once.
    ///
    /// Setting this to a larger value will allow more atoms to be formatted at once, at the cost
    /// of higher memory consumption. Setting this to a smaller value will lower the memory
    /// footprint at a possible cost of efficiency.
    ///
    /// Note that depending on the number of residues in the template box, the provided value may
    /// be overshot by at most that number.
    #[arg(long, default_value_t = 10000000)]
    pub buffer_size: usize,
}

#[derive(Debug, Default, Clone, Copy, ValueEnum, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub enum WaterType {
    #[default]
    Martini,
    Tip3P,
}

impl WaterType {
    pub const fn default_cutoff(&self) -> f32 {
        match self {
            WaterType::Martini => 0.43,
            WaterType::Tip3P => 0.28,
        }
    }

    pub const fn default_solvent_cutoff(&self) -> f32 {
        match self {
            WaterType::Martini => 0.21,
            WaterType::Tip3P => 0.21,
        }
    }
}

#[derive(Debug, Default, Clone, ValueEnum, PartialEq, Eq)]
pub enum BoundaryMode {
    /// Cut at the structure box size and remove solvent residues that overlap with the periodic
    /// neighbors.
    #[default]
    Cut,
    /// If necessary, grow the box of the output structure to fit a whole number of template boxes.
    Grow,
}

#[derive(Debug, Default, Clone, ValueEnum, PartialEq, Eq)]
pub enum PeriodicMode {
    /// Include the periodic images of the structure when checking whether a solvent spot is
    /// occupied.
    #[default]
    Periodic,
    /// Ignore the periodic images of the structure for the solvent placement check.
    Ignore,
    /// Treat any structure atoms outside of the output box as an error and exit.
    Deny,
}

#[derive(Debug, Clone)]
pub struct SubstituteConfig {
    name: String,
    quantity: Quantity,
}

impl SubstituteConfig {
    /// Bake into a [`Substitute`] according to a final volume in nm³ and some number of valid
    /// solvent beads.
    pub fn bake(self, volume: f64, solvent_beads: u64) -> Substitute {
        Substitute {
            name: self.name,
            number: self.quantity.number(volume, solvent_beads),
        }
    }
}

impl FromStr for SubstituteConfig {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split(':');
        let name = words
            .next()
            .ok_or("expected a name, then a colon, followed by a quantity".to_string())?
            .to_string();
        let quantity = words
            .next()
            .ok_or("expected a quantifier after the colon".to_string())?
            .parse()?;
        Ok(Self { name, quantity })
    }
}

#[derive(Debug, Clone)]
pub enum Charge {
    Known(i64),
    Auto,
}

#[derive(Debug, Clone)]
pub struct ChargeConfig {
    charge: Charge,
    positive: String,
    negative: String,
}

impl FromStr for ChargeConfig {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split(':');
        let charge = words.next().ok_or("expected a charge")?;
        let charge = match charge {
            "auto" => Charge::Auto,
            n => Charge::Known(n.parse::<i64>().map_err(|err| err.to_string())?),
        };
        let (positive, negative) = if let Some(names) = words.next() {
            names
                .split_once(',')
                .ok_or("expected two substitute names separated by a comma")?
        } else {
            ("NA", "CL")
        };
        Ok(Self {
            charge,
            positive: positive.to_string(),
            negative: negative.to_string(),
        })
    }
}

impl ChargeConfig {
    /// Bake the [`ChargeConfig`] into a [`Substitute`] if applicable.
    ///
    /// If the charge is 0, no `Substitute` is returned since there is no charge to neutralize.
    pub fn bake<P: AsRef<std::path::Path> + std::fmt::Debug>(
        self,
        topol_path: Option<P>,
    ) -> Result<Option<Substitute>> {
        // Choose the name of the ion to neutralize with. If the charge to compensate is negative,
        // we compensate with the positive ion substitute, and vice versa.
        let charge = match self.charge {
            Charge::Known(known) => known,
            Charge::Auto => {
                let topol_path =
                    topol_path.context("automatic charge determination requires a topology")?;
                let start = std::time::Instant::now();
                let charge = determine_system_charge(&topol_path).context(format!(
                    "could not determine system charge from topology file at {topol_path:?}"
                ))?;
                let duration = start.elapsed().as_secs_f32();
                eprintln!("According to the topology, the total system charge is {charge}.");
                eprintln!("Determining charges took {duration:.3} s.");
                charge
            }
        };
        let name = match charge {
            ..0 => self.positive,
            0 => return Ok(None),
            1.. => self.negative,
        };

        Ok(Some(Substitute {
            name,
            number: charge.unsigned_abs(),
        }))
    }
}

fn determine_system_charge<P: AsRef<std::path::Path> + std::fmt::Debug>(
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

    // "Using [ molecules ] without having used [ system ] before is meaningless and generates a
    // warning."
    skip_to_section(&mut lines, "system");
    // The next line is the system title.
    let (_ln, _title) = lines.next().context("expected system title")?;

    // "After [ system ] the only allowed directive is [ molecules ]."
    skip_to_section(&mut lines, "molecules");
    // The next lines contain the list of molecules in the system.
    for (ln, line) in lines {
        // Parse a molecule.
        let Some((name, number)) = line.split_once(char::is_whitespace) else {
            bail!("expected molecule declaration at line {ln}, found {line:?}");
        };
        let number: MoleculeNumber = number
            .trim_start()
            .parse()
            .with_context(|| format!("expected an integer, found {number:?} at line {ln}"))?;
        *molecules.entry(name).or_default() += number;
    }

    // Now we go through the top file a second time, dilligently following all includes and
    // processing the relevant [ molecule_type ] sections.
    // That means that we are going to treat the root top file as if it were any normal itp file,
    // and we follow every include in a depth-first manner.

    // // Create a map to the mutable charges that is itself immutable, such that its contents are not
    // // suddenly rearranged.
    // let relevant_molecules: Box<[MoleculeName]> = molecules.keys().copied().collect();
    // let mut moleculetypes: Box<[_]> = relevant_molecules
    //     .into_iter()
    //     .map(|name| (name, Default::default()))
    //     .collect();

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
        line.strip_prefix('[')
            .and_then(|rest| rest.strip_suffix(']'))
            .map(str::trim)
    }

    fn parse_path_str<'s>(s: &'s str, ln: usize, line: &'s str) -> Result<&'s str> {
        let front = s.trim_start().strip_prefix('"').with_context(|| {
            format!("expected the start of a quoted file path at line {ln}, found {line:?}")
        })?;
        front.strip_suffix('"').with_context(|| {
            format!("expected the end of a quoted file path at line {ln}, found {line:?}")
        })
    }

    struct Molecule {
        name: String,
        charge: f64,
    }

    impl Molecule {
        fn new(name: String) -> Self {
            Self {
                name,
                charge: Default::default(),
            }
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
                                .with_context(|| format!("could not parse {full_path:?}"))?;
                        }

                        // Change the conditionals stack.
                        "ifdef" => {
                            let symbol = rest
                                .split_whitespace()
                                .next()
                                .context("expected a definition name")?;
                            // We go a level deeper.
                            let current = self.defined.contains(symbol);
                            let parent_active = self.is_active();
                            self.stack.push(current && parent_active);
                        }
                        "ifndef" => {
                            let symbol = rest
                                .split_whitespace()
                                .next()
                                .context("expected a definition name")?;
                            // We go a level deeper.
                            let current = !self.defined.contains(symbol);
                            let parent_active = self.is_active();
                            self.stack.push(current && parent_active);
                        }
                        "else" => {
                            let prev = self.stack.pop().context(
                                "encountered #else directive without preceding #ifdef or #ifndef",
                            )?;
                            // Invert our activity state in the else branch.
                            let parent_active = self.is_active();
                            self.stack.push(parent_active && !prev);
                        }
                        "endif" => {
                            self.stack.pop().context("popped too far!")?; // Pop down a level.
                        }

                        // Important to not consider any #define or #undef directives when we are
                        // not active.
                        "define" if self.is_active() => {
                            let symbol = rest
                                .split_whitespace()
                                .next()
                                .context("expected a definition name")?;
                            self.defined.insert(symbol.to_string());
                        }
                        "undef" if self.is_active() => {
                            let symbol = rest
                                .split_whitespace()
                                .next()
                                .context("expected a definition name")?;
                            self.defined.insert(symbol.to_string());
                        }

                        unknown => {
                            eprintln!(
                                "UNKNOWN DIRECTIVE: #{unknown} at {ln} in {path:?}:\n{line:?}"
                            );
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
                        format!(
                            "could not parse atom parameters line at {ln} in {path:?}: {line:?}"
                        )
                    })?,
                    State::MoleculeType => {
                        if let Some(name) = line.split_whitespace().next() {
                            self.current = Some(Molecule::new(name.to_string()));
                        }
                    }
                    State::Atoms => self.parse_molecule_atom(line).with_context(|| {
                        let name = self
                            .current
                            .as_ref()
                            .map(|curr| curr.name.as_str())
                            .unwrap_or("?");
                        format!(
                            "could not parse atom line at {ln} in {path:?} \
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
                bail!("too few fields in atom type definition");
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
                bail!("the particle type for this atom line cannot be reliably determined")
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
                bail!("cannot read an atoms section without a prior [ moleculetype ] directive")
            };

            // Fields: nr, atom type, residue number, residue name, atom name, charge group number, q(e), m(u).
            let atom_type_name = fields.nth(1).context("expected atom type name")?;
            // If the moleculetype atom contains a charge, that overrides the charge set in the
            // forcefield parameters for that atomtype.
            let charge = if let Some(charge) = fields.nth(4) {
                charge
                    .parse()
                    .with_context(|| format!("could not parse charge value in atom line"))?
            } else if let Some(&charge) = self.atoms.get(atom_type_name) {
                charge
            } else {
                // I think this would be considered an error in the eyes of GROMACS as well.
                bail!(
                    "atom type {atom_type_name} was not previously defined \
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
    }

    let moleculetypes = {
        let mut context = Context::default();
        context
            .parse(topol_path.as_ref())
            .with_context(|| format!("problem while parsing main topology file {topol_path:?}"))?;
        context.wrap_up();
        context.molecules
    };

    // Finally, the accounting. We now know how many of each moleculetype we have, and what the
    // charge of each individual instance of that moleculetype is. Now we compute the total charge.
    let mut total_charge = 0.0;
    for (name, number) in molecules {
        let Some(charge): Option<f64> = moleculetypes
            .iter()
            .find_map(|m| if m.name == name { Some(m.charge) } else { None })
        else {
            bail!("no moleculetype was defined for {name:?}");
        };

        // Report any fractional charges of individual structures as warnings. This makes it much
        // easier to identify issues since we have a hard error on fractional charges below.
        if f64::fract(charge.round() - charge).abs() > FRACTIONAL_CHARGE_THRESHOLD {
            eprintln!("WARNING: Molecule {name:?} has a fractional charge of {charge}.");
        }

        total_charge += number as f64 * charge;
    }

    // TODO: This threshold can be set with greater care maybe. Just a guess.
    const FRACTIONAL_CHARGE_THRESHOLD: f64 = 0.01;
    if f64::fract(total_charge.round() - total_charge).abs() > FRACTIONAL_CHARGE_THRESHOLD {
        bail!("the system charge should be an integer value, but it is {total_charge}");
    }
    let total_charge = total_charge.round() as i64;

    Ok(total_charge)
}

#[derive(Debug, Clone)]
enum Quantity {
    /// A fixed number.
    Number(u64),
    /// A fraction of the valid solvent beads.
    Ratio(f64),
    /// Molarity in mol/L.
    Molarity(f64),
}

impl FromStr for Quantity {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(molarity) = s.strip_suffix('M') {
            // Molarity.
            return Ok(Self::Molarity(
                molarity.parse::<f64>().map_err(|err| err.to_string())?,
            ));
        }

        if let Ok(number) = s.parse::<u64>() {
            return Ok(Self::Number(number));
        }

        if let Ok(ratio) = s.parse::<f64>() {
            if !(0.0..=1.0).contains(&ratio) {
                return Err(format!(
                    "ratio value must satisfy `1.0 >= r >= 0.0`, found {ratio}"
                ));
            }
            return Ok(Self::Ratio(ratio));
        }

        Err(format!("could not parse quantity {s:?}"))
    }
}

impl Quantity {
    /// Given some volume in nm³ and some number of initial valid solvent beads, determine the
    /// number of beads according to this [`Quantity`].
    fn number(self, volume: f64, solvent_beads: u64) -> u64 {
        const N_AVOGADRO: f64 = 6.0221415e23; // per mol.
        let volume_liter = volume * 1e-24;
        let number = match self {
            Quantity::Number(number) => number,
            Quantity::Ratio(ratio) => (solvent_beads as f64 * ratio).round() as u64,
            Quantity::Molarity(molarity) => (volume_liter * molarity * N_AVOGADRO).round() as u64,
        };
        assert!(
            number <= solvent_beads,
            "the number of substitutions that are specified ({number}) exceeds the number of \
                solvent positions that can be substituted ({solvent_beads})"
        );
        number
    }
}

pub struct Substitute {
    pub name: String,
    pub number: u64,
}

#[derive(ValueEnum, Default, Clone, Debug)]
pub enum SortBehavior {
    #[default]
    Size,
    RevSize,
    Alphabetical,
    RevAlphabetical,
    No,
}
