use std::borrow::Cow;
use std::collections::HashSet;
use std::io::{BufWriter, Read, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand, command};

use bentopy::core::config::{Config, Segment, legacy};
use bentopy::core::version::VERSION;

const BIN_NAME: &str = env!("CARGO_BIN_NAME");

#[derive(Debug, Parser)]
#[command(about, version = VERSION)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Example {
        #[arg(short, long, default_value = "example.bent")]
        output: PathBuf,
    },
    // Create {
    //     #[arg(short, long)]
    //     output: PathBuf,
    // },
    Check {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value_t)]
        verbose: bool,
    },
    Convert {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
    },
}

fn main() -> Result<()> {
    let Args { command } = Args::parse();
    match command {
        Command::Example { output } => example(output),
        // Command::Create { output } => todo!(),
        Command::Check { input, verbose } => check(input, verbose),
        Command::Convert { input, output } => convert(input, output),
    }
}

fn example(output: PathBuf) -> Result<()> {
    let example = include_str!("example.bent");
    let mut file = std::fs::File::create(&output)?;
    writeln!(file, "# Created by {BIN_NAME}, version {VERSION}.")?;
    file.write_all(example.as_bytes())?;
    eprintln!("Wrote example file to {output:?}.");
    Ok(())
}

fn check(input: PathBuf, verbose: bool) -> Result<()> {
    let mut file = std::fs::File::open(&input)?;
    let mut s = String::new();
    file.read_to_string(&mut s)?;

    // Try to parse the config.
    let config = Config::parse_bent(&input.to_string_lossy(), &s)
        .context(format!("could not parse {input:?} as a bentopy input file"))?;
    eprintln!("Successfully parsed {input:?}.");

    // If desired, print a crazy expanded debug struct.
    if verbose {
        eprintln!("Printing debug representation of this configuration:");
        println!("{config:#?}");
    }

    // TODO: Write out a nice little report explaining the system. Maybe this is actually a
    // bentopy-init explain thing but whatever. It's a good idea.

    // TODO: Implement these.
    // eprintln!("NOTE: This is a work-in-progress utility. Notable lacking lints are:");
    // eprintln!("      - check if the referenced structures exist,");
    // eprintln!("      - check if the compartment combinations are well-formed,");
    // eprintln!("      - check if the constraint combinations are well-formed,");
    // eprintln!("      - check if the referenced masks actually exist,");
    // eprintln!("      - check if the referenced masks are of a congruent size,");
    // eprintln!("      - check if the referenced itp includes actually exist,");
    // eprintln!("      - check if the segment names exist in the referenced itps,");
    // eprintln!("      - check if the analytical compartment geometries lie within the dimensions,");
    // eprintln!("      - check for orphan rules and compartments,");
    // eprintln!("      - check if the rules apply within the dimensions.");

    // Now, we actually check if it is semantically correct.
    let mut problems: usize = 0;
    // There should be at least one compartment.
    if config.compartments.is_empty() {
        println!("Error: No compartments declared. At least one compartment should be declared.");
        problems += 1;
    }
    // There should be no duplicate compartment names.
    let mut compartment_ids = HashSet::with_capacity(config.compartments.len());
    for compartment in &config.compartments {
        let id = compartment.id.as_str();
        if !compartment_ids.insert(id) {
            println!("Warning: Compartment id {id} is used multiple times.");
            problems += 1;
        }
    }
    // There should be at least one segment.
    if config.compartments.is_empty() {
        println!("Warning: No segments declared. At least one segment should be declared.");
        problems += 1;
    }
    // Segments should refer to compartments that exist.
    for segment in &config.segments {
        let name = segment.name();
        for id in &segment.compartment_ids {
            if !compartment_ids.contains(id.as_str()) {
                println!("Error: Segment {name} refers to an undeclared compartment: {id}.");
                problems += 1;
            }
        }
    }
    // There should be no duplicate rule names.
    let mut rules = HashSet::with_capacity(config.constraints.len());
    for constraint in &config.constraints {
        let id = constraint.id.as_str();
        if !rules.insert(id) {
            println!("Warning: Constraint name {id} is used multiple times.");
            problems += 1;
        }
    }
    // Segments should refer to rules that exist.
    for segment in &config.segments {
        let name = segment.name();
        for id in segment.rules.iter() {
            if !rules.contains(id.as_str()) {
                println!("Error: Segment {name} refers to an undeclared constraint: {id}.");
                problems += 1;
            }
        }
    }

    // Report the number of problems.
    match problems {
        0 => println!("It appears there are no problems with this file."),
        1 => println!("Detected one problem."),
        n => println!("Detected {n} problems."),
    }

    Ok(())
}

fn convert(input: PathBuf, output: PathBuf) -> Result<()> {
    enum Kind {
        Bent,
        Json,
    }

    let kind = match input.extension().and_then(|s| s.to_str()) {
        Some("bent") => Kind::Bent,
        Some("json") => Kind::Json,
        _ => bail!(
            "Cannot read {input:?}. \
                The convert subcommand only supports reading .bent and .json files.",
        ),
    };

    if output.extension().is_some_and(|s| s != "bent") {
        bail!(
            "Cannot write configuration to {output:?}. \
                The convert subcommand only supports writing .bent files.",
        );
    }

    let mut file = std::fs::File::open(&input)?;
    let mut s = String::new();
    file.read_to_string(&mut s)?;

    let path = &input.to_string_lossy();
    let config = match kind {
        Kind::Bent => {
            let config = Config::parse_bent(path, &s)
                .context(format!("could not parse {path:?} as a bentopy input file"))?;
            eprintln!("Successfully parsed {path:?}.");
            config
        }
        Kind::Json => {
            let config: legacy::Config = serde_json::from_str(&s).context(format!(
                "could not parse {path:?} as a legacy json input file"
            ))?;
            eprintln!("Successfully parsed {path:?} (legacy input file).");
            eprint!("Attempting to convert to new input configuration format... ");
            let config = config.into();
            eprintln!("Done.");
            config
        }
    };

    eprintln!("Writing to {output:?}.");
    let out = std::fs::File::create(&output)?;
    let mut out = BufWriter::new(out);
    writeln!(out, "# Converted {input:?} to {output:?} using {BIN_NAME}.")?;
    bentopy::core::config::bent::write(&config, &mut out)?;

    Ok(())
}

// TODO: Make this into a method on Segment itself. It would also be nice in bent::writer.
trait Name {
    fn name(&self) -> Cow<'_, str>;
}

impl Name for Segment {
    fn name(&self) -> Cow<'_, str> {
        let name = &self.name;
        match &self.tag {
            Some(tag) => Cow::Owned(format!("{name}:{tag}")),
            None => Cow::Borrowed(name),
        }
    }
}
