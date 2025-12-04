use std::io::{BufWriter, Read, Write};
use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand, command};

use bentopy::core::config::Config;

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
    let mut file = std::fs::File::open(&input)?;
    let mut s = String::new();
    file.read_to_string(&mut s)?;

    let config = Config::parse_bent(&input.to_string_lossy(), &s)?;
    eprintln!("Successfully parsed {input:?}.");
    eprintln!("Writing to {output:?}.");
    let out = std::fs::File::create(output)?;
    let mut out = BufWriter::new(out);
    bentopy::core::config::bent::write(&config, &mut out)?;

    Ok(())
}
