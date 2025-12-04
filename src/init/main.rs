use std::io::{BufWriter, Read, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand, command};

use bentopy::core::config::{Config, legacy};

#[derive(Debug, Parser)]
#[command(about, version = bentopy::core::version::VERSION)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Example {
        #[arg(short, long)]
        output: PathBuf,
    },
    // Create {
    //     #[arg(short, long)]
    //     output: PathBuf,
    // },
    Check {
        #[arg(short, long)]
        input: PathBuf,
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
        Command::Check { input } => check(input),
        Command::Convert { input, output } => convert(input, output),
    }
}

fn example(output: PathBuf) -> Result<()> {
    let example = include_str!("example.bent");
    let mut file = std::fs::File::create(&output)?;
    file.write_all(example.as_bytes())?;
    Ok(())
}

fn check(input: PathBuf) -> Result<()> {
    let mut file = std::fs::File::open(&input)?;
    let mut s = String::new();
    file.read_to_string(&mut s)?;

    let config = Config::parse_bent(&input.to_string_lossy(), &s)?;
    eprintln!("Successfully parsed {input:?}.");
    println!("{config:#?}");

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
            eprintln!("Attempting to convert to new input configuration format...");
            let config = config.into();
            eprintln!("Done!");
            config
        }
    };

    eprintln!("Writing to {output:?}.");
    let out = std::fs::File::create(output)?;
    let mut out = BufWriter::new(out);
    bentopy::core::config::bent::write(&config, &mut out)?;

    Ok(())
}
