use anyhow::{Context, Result, bail};

use crate::core::config::{Center, Compartment, Config, Mask, Point, Quantity, Segment, Shape};

pub fn parse_config(s: &str) -> Result<Config> {
    let mut config = Config {
        title: "System".to_owned(),
        seed: None,
        bead_radius: 0.20,
        max_tries_mult: 1000,
        max_tries_rot_div: 100,
        dimensions: [10.0; 3],
        resolution: 0.5,
        periodic: true,
        compartments: Default::default(),
        topol_includes: Default::default(),
        segments: Default::default(),
    };

    let mut lines = s.lines().enumerate().peekable();
    while let Some((ln, line)) = lines.next() {
        let Some(line) = strip_comments(line) else {
            continue;
        };

        // At this point, any remaining line has no surrounding spaces nor trailing comments.
        if let Some(potential_header) = line.strip_prefix('[')
            && let Some(header) = potential_header.strip_suffix(']')
        {
            // A header is surrounded by brackets.
            let header = header.trim(); // "Tighten up those lines!"
            let declaration_parser = match header {
                "general" => parse_general_declaration,
                "space" => parse_space_declaration,
                "compartments" => parse_compartments_declaration,
                "includes" => parse_includes_declaration,
                "segments" => parse_segments_declaration,
                unknown => {
                    bail!("encountered an unknown config header on line {ln}: {unknown:?}")
                }
            };
            parse_section(&mut lines, &mut config, declaration_parser)
                .context(format!("could not parse the {header} section"))?;
        } else {
            // Otherwise, we're dealing with an orphan line.
            bail!("encountered a declaration not under a header at line {ln}: {line:?}")
        }
    }

    Ok(config)
}

/// Strip any comments.
///
/// Returns [`Some`] line if the line is not empty. If the line is empty,
/// this function returns [`None`].
fn strip_comments(line: &str) -> Option<&str> {
    // Strip any comments.
    let line = match line.split_once('#') {
        Some((line, _comment)) => line,
        None => line,
    }
    .trim();
    if line.is_empty() {
        // Skip empty lines and line comments.
        return None;
    }
    Some(line)
}

fn parse_vec3(s: &str) -> Result<Point> {
    let mut values = s.split(',');
    let mut parse_next = |name| {
        values
            .next()
            .context("expected a floating point number")?
            .parse()
            .with_context(|| format!("expected {name} value as a floating point number"))
    };
    let x = parse_next("x")?;
    let y = parse_next("y")?;
    let z = parse_next("z")?;

    Ok([x, y, z])
}

fn parse_section<'a, I>(
    lines: &mut std::iter::Peekable<I>,
    config: &mut Config,
    declaration_parser: fn(&mut Config, usize, &str) -> Result<()>,
) -> Result<()>
where
    I: Iterator<Item = (usize, &'a str)>,
{
    loop {
        // First, we check if we are running into the next header or the end of the file.
        // We leave that to be handled after we return.
        match lines.peek() {
            // Encountered a header. Exiting.
            Some((_ln, line)) if line.trim_start().starts_with('[') => break,
            // We are at the end. Exiting.
            None => break,
            _ => {}
        }

        // Let's take the next line now.
        let (ln, line) = lines.next().unwrap(); // We know it exists.
        let Some(line) = strip_comments(line) else {
            continue;
        };

        // Now we know that we are dealing with a declaration line.
        declaration_parser(config, ln, line)
            .context(format!("could not parse declaration {line:?} on line {ln}"))?
    }

    Ok(())
}

fn parse_general_declaration<'a>(
    config: &mut Config,
    ln: usize,
    line: &str,
) -> std::result::Result<(), anyhow::Error> {
    let (keyword, value) = parse_keyword_value(line, ln, "<keyword> <value>")?;
    match keyword {
        "title" => config.title = value.to_owned(),
        "seed" => {
            config.seed = Some(
                value
                    .parse()
                    .context(format!("could not parse seed with value {value:?}"))?,
            )
        }
        "bead-radius" => {
            config.bead_radius = value
                .parse()
                .context(format!("could not parse radius with value {value:?}"))?
        }
        "max-tries-mult" => {
            config.max_tries_mult = value.parse().context(format!(
                "could not parse max_tries_mult with value {value:?}"
            ))?
        }
        "max-tries-rot-div" => {
            config.max_tries_rot_div = value.parse().context(format!(
                "could not parse max_tries_rot_div with value {value:?}"
            ))?
        }
        keyword => bail!("unknown keyword {keyword:?} in 'general' section on line {ln}"),
    }

    Ok(())
}

/// Parse a line into a whitespace-separated keyword-value pair.
fn parse_keyword_value<'s>(
    line: &'s str,
    ln: usize,
    exp: &'static str,
) -> Result<(&'s str, &'s str)> {
    let Some((keyword, value)) = line.split_once(char::is_whitespace) else {
        bail!("expected a declaration of the form '{exp}' on line {ln}, but found {line:?}")
    };
    Ok((keyword, value.trim_start()))
}

fn parse_space_declaration<'a>(
    config: &mut Config,
    ln: usize,
    line: &str,
) -> std::result::Result<(), anyhow::Error> {
    let (keyword, value) = parse_keyword_value(line, ln, "<keyword> <value>")?;
    match keyword {
        "dimensions" => {
            config.dimensions = parse_vec3(value).context("could not parse dimensions")?
        }
        "resolution" => config.resolution = value.parse().context("could not parse resolution")?,
        "periodic" => config.periodic = value.parse().context("could not parse periodic (bool)")?,
        keyword => bail!("unknown keyword {keyword:?} in 'space' section on line {ln}"),
    }

    Ok(())
}

fn parse_compartments_declaration<'a>(
    config: &mut Config,
    ln: usize,
    line: &str,
) -> std::result::Result<(), anyhow::Error> {
    let (id, mask_declaration) = parse_keyword_value(line, ln, "<id> <mask-declaration>")?;
    let compartment = Compartment {
        id: id.to_owned(),
        mask: parse_mask(mask_declaration).context(format!(
            "could not parse mask declaration {mask_declaration:?} on line {ln}"
        ))?,
    };
    config.compartments.push(compartment);

    Ok(())
}

fn parse_mask(d: &str) -> Result<Mask> {
    let Some((keyword, rest)) = d.split_once(char::is_whitespace) else {
        bail!("expected mask declaration")
    };
    let mask = match keyword {
        "from" => Mask::Voxels(rest.into()),
        "as" => Mask::Shape(
            parse_shape(rest).context(format!("could not parse shape declaration {rest:?}"))?,
        ),
        unknown => bail!("expected 'as <shape>' or 'from <path>', found {unknown:?}"),
    };

    Ok(mask)
}

fn parse_shape(d: &str) -> Result<Shape> {
    let Some((keyword, rest)) = d.split_once(char::is_whitespace) else {
        bail!("expected shape declaration")
    };
    let shape = match keyword {
        "sphere" => {
            let Some(rest) = rest.strip_prefix("at") else {
                bail!("expected 'at <center> with <radius-or-diameter> <size>'")
            };
            let mut words = rest.split_whitespace();

            let center = match words
                .next()
                .context("expected a center definition (point or 'center')")?
            {
                "center" => Center::Center,
                point => Center::Point(parse_vec3(point)?),
            };

            let Some("with") = words.next() else {
                bail!(
                    "expected sphere radius definition 'with radius <size>' or 'with diameter <size>'"
                );
            };
            let kind = words.next().context("expected 'radius' or 'diameter'")?;
            let size = words
                .next()
                .context("expected size as floating point number")?
                .parse()
                .context("could not parse size")?;
            let radius = match kind {
                "radius" => size,
                "diameter" => size * 0.5,
                unknown => {
                    bail!(
                        "unknown sphere size keyword {unknown:?}, expected 'radius' or 'diameter'"
                    )
                }
            };
            Shape::Sphere { center, radius }
        }
        "cuboid" => {
            let Some(rest) = rest.strip_prefix("from") else {
                bail!("expected 'from <start> to <end>'")
            };
            let mut words = rest.split_whitespace();
            let start = parse_vec3(
                words
                    .next()
                    .context("expected a cuboid start point definition")?,
            )?;
            let Some("to") = words.next() else {
                bail!("expected 'to <end>'")
            };
            let end = parse_vec3(
                words
                    .next()
                    .context("expected a cuboid end point definition")?,
            )?;
            Shape::Cuboid { start, end }
        }
        "combination" => Shape::Combination {
            expression: rest
                .parse()
                .context("could not parse compartment combination expression")?,
        },
        unknown => bail!("unknown shape keyword {unknown:?}, expected 'sphere' or 'cuboid'"),
    };

    Ok(shape)
}

fn parse_includes_declaration<'a>(
    config: &mut Config,
    _ln: usize,
    path: &str,
) -> std::result::Result<(), anyhow::Error> {
    // TODO: Do proper string quoting and escaping!
    config.topol_includes.push(path.into());
    Ok(())
}

fn parse_segments_declaration<'a>(
    config: &mut Config,
    ln: usize,
    line: &str,
) -> std::result::Result<(), anyhow::Error> {
    let (name, segment_declaration) =
        parse_keyword_value(line, ln, "<name> <segment-declaration>")?;
    let (tag, rest) = if let Some(tag) = segment_declaration.strip_prefix('(') {
        let Some((tag, rest)) = tag.split_once(')') else {
            bail!("expected closing parenthesis after tag");
        };
        (Some(tag.to_owned()), rest)
    } else {
        (None, segment_declaration)
    };

    let mut words = rest.split_whitespace();
    let Some("from") = words.next() else {
        bail!("expected 'from'");
    };
    let Some(path) = words.next() else {
        bail!("expected a structure path");
    };
    let Some("in") = words.next() else {
        bail!("expected 'in'");
    };
    let Some(compartment_ids) = words.next() else {
        bail!("expected compartment names");
    };
    let Some("at") = words.next() else {
        bail!("expected 'at'");
    };
    let Some(quantity) = words.next() else {
        bail!("expected quantity (copy number or molarity)");
    };

    let segment = Segment {
        name: name.to_owned(),
        tag,
        path: path.into(),
        compartment_ids: compartment_ids
            .split(',')
            .map(ToString::to_string)
            .collect(),
        quantity: parse_quantity(quantity).context("error while parsing quantity")?,
    };
    config.segments.push(segment);

    Ok(())
}

fn parse_quantity(s: &str) -> Result<Quantity> {
    let quantity = if let Some(molarity) = s.strip_suffix("M") {
        let molarity = molarity
            .parse()
            .context(format!("could not parse {s:?} as molarity"))?;
        Quantity::Concentration(molarity)
    } else {
        let number = s
            .parse()
            .context(format!("could not parse {s:?} as copy number"))?;
        Quantity::Number(number)
    };
    Ok(quantity)
}
