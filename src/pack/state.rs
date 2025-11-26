use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{Context, bail};
pub(crate) use glam::{EulerRot, Mat3, Vec3};
use rand::{RngCore, SeedableRng};

pub use bentopy::core::config::CompartmentID;
use bentopy::core::config::legacy::{
    Compartment as LegacyConfigCompartment, Config as LegacyConfig, Mask as LegacyConfigMask,
    RuleExpression, TopolIncludes,
};

use crate::args::{Args, RearrangeMethod};
use crate::mask::Mask;
use crate::rules::{self, ParseRuleError, Rule};
pub use crate::state::compartment::Compartment;
use crate::state::segment::Segment;
pub use crate::state::space::Space;
use crate::structure::load_molecule;

mod combinations;
mod compartment;
mod mask;
mod pack;
mod segment;
mod space;

const ORDER: EulerRot = EulerRot::XYZ;

pub type Size = [f32; 3];
pub type Rotation = Mat3;
pub type Voxels = Mask;
pub type Rng = rand::rngs::StdRng; // TODO: Is this the fastest out there?

pub struct Output {
    pub title: String,
    pub path: PathBuf,
    pub topol_includes: Option<TopolIncludes>,
}

pub struct General {
    pub seed: u64,
    pub max_tries_multiplier: u64,
    pub max_tries_per_rotation_divisor: u64,
    pub bead_radius: f32,
}

pub struct State {
    pub general: General,
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,

    pub rng: Rng,
    pub verbose: bool,
    pub summary: bool,
}

impl State {
    /// Set up the [`State`] given the [command line arguments](Args) and input file
    /// [configuration](LegacyConfig).
    pub fn new(args: Args, config: LegacyConfig) -> anyhow::Result<Self> {
        // Read values from the general section of the config. If a command line argument is given,
        // it overwrites the config value. (And the deprecated env vars have the highest priority.)

        // If no seed is provided, use a random seed.
        let seed = args
            .seed
            .or(config.general.seed)
            .unwrap_or_else(|| Rng::from_os_rng().next_u64());
        let rng = Rng::seed_from_u64(seed);

        let bead_radius = args.bead_radius.unwrap_or(config.general.bead_radius);

        // Determine the max_tries parameters.
        let max_tries_multiplier = if let Ok(s) = std::env::var("BENTOPY_TRIES") {
            let n = s.parse().with_context(|| {
                format!("Max tries multiplier should be a valid unsigned integer, found {s:?}")
            })?;
            eprintln!("\tMax tries multiplier set to {n}.");
            eprintln!(
                "\tWARNING: Setting max_tries_mult using the BENTOPY_TRIES environment variable will be deprecated."
            );
            eprintln!("\t         Use --max-tries-mult instead.");
            n
        } else {
            args.max_tries_mult.unwrap_or(config.general.max_tries_mult)
        };

        let max_tries_per_rotation_divisor = if let Ok(s) = std::env::var("BENTOPY_ROT_DIV") {
            let n = s.parse().with_context(|| {
                format!("Rotation divisor should be a valid unsigned integer, found {s:?}")
            })?;
            eprintln!("\tMax tries per rotation divisor set to {n}.");
            eprintln!(
                "\tWARNING: Setting max_tries_divisor using the BENTOPY_ROT_DIV environment variable will be deprecated."
            );
            eprintln!("\t         Use --max-tries-rot-div instead.");
            n
        } else {
            args.max_tries_rot_div
                .unwrap_or(config.general.max_tries_rot_div)
        };

        let verbose = args.verbose;

        // Space.
        let dimensions = config
            .space
            .size
            .map(|d| (d / config.space.resolution) as u64);
        let resolution = config.space.resolution;
        eprintln!("Setting up compartments...");
        let (predefined, combinations): (Vec<_>, Vec<_>) = config
            .space
            .compartments
            .into_iter()
            .partition(LegacyConfigCompartment::is_predefined);
        let mut compartments: Vec<Compartment> = predefined
            .into_iter()
            .map(|comp| -> anyhow::Result<_> {
                let mask = match comp.mask {
                    LegacyConfigMask::Shape(shape) => {
                        if verbose {
                            eprintln!("\tConstructing a {shape} mask...");
                        }
                        Mask::legacy_create_from_shape(shape, dimensions, None, None)
                    }
                    LegacyConfigMask::Analytical {
                        shape,
                        center,
                        radius,
                    } => {
                        if verbose {
                            eprintln!("\tConstructing a {shape} mask...");
                        }
                        let center = center.map(|c| (Vec3::from_array(c) / resolution).as_uvec3());
                        let radius = radius.map(|r| (r / resolution) as u32);
                        Mask::legacy_create_from_shape(shape, dimensions, center, radius)
                    }
                    LegacyConfigMask::Voxels { path } => {
                        if verbose {
                            eprintln!("\tLoading mask from {path:?}...");
                        }
                        Mask::load_from_path(&path)
                            .with_context(|| format!("Failed to load mask {path:?}"))?
                    }
                    LegacyConfigMask::Combination(_) => {
                        unreachable!() // We partitioned the list above.
                    }
                };
                let c = Compartment {
                    id: comp.id,
                    mask,
                    distance_masks: Default::default(),
                };
                Ok(c)
            })
            .collect::<anyhow::Result<_>>()?;
        for combination in combinations {
            let LegacyConfigMask::Combination(ces) = combination.mask else {
                unreachable!() // We partitioned the list above.
            };
            let baked = Compartment {
                id: combination.id,
                mask: combinations::execute(&ces, &compartments)?,
                distance_masks: Default::default(),
            };
            compartments.push(baked);
        }
        let space = Space {
            size: config.space.size,
            dimensions,
            resolution,
            compartments,
            periodic: config.space.periodic,

            global_background: Mask::new(dimensions),
            session_background: Mask::new(dimensions),
            previous_compartments: None,
            previous_rules: None,
        };

        // Segments.
        eprintln!("Loading segment structures...");
        let segments = {
            let mut segments: Vec<_> = config
                .segments
                .into_iter()
                .map(|seg| -> Result<_, _> {
                    if verbose {
                        eprintln!("\tLoading {:?}...", &seg.path);
                    }
                    let name = seg.name;
                    let tag = seg.tag;
                    match tag.as_ref().map(String::len) {
                        Some(0) => eprintln!("WARNING: The tag for segment '{name}' is empty."),
                        Some(6.. ) => eprintln!("WARNING: The tag for segment '{name}' is longer than 5 characters, and may be truncated when the placement list is rendered."),
                        _ => {} // Nothing to warn about.
                    }
                    let path = seg.path;
                    let structure = load_molecule(&path).with_context(|| format!("Failed to open the structure file for segment '{name}' at {path:?}"))?;
                    let rules = seg
                        .rules
                        .iter()
                        .map(parse_rule)
                        .collect::<Result<Vec<Rule>, _>>()?;
                    let initial_rotation = {
                        let [ax, ay, az] = seg.initial_rotation.map(f32::to_radians);
                        Rotation::from_euler(ORDER, ax, ay, az)
                    };
                    Ok(Segment {
                        name,
                        tag,
                        quantity: seg.quantity,
                        compartments: seg.compartments,
                        rules,
                        path,
                        rotation_axes: seg.rotation_axes,
                        structure,
                        initial_rotation,
                        rotation: Rotation::IDENTITY,
                        voxels: None,
                    })
                })
                .collect::<anyhow::Result<_>>()?;

            let method = args.rearrange;
            if let RearrangeMethod::None = method {
                eprint!("Segments were not rearranged.");
            } else {
                eprint!("Rearranging segments according to the {method:?} method... ");
                match method {
                    RearrangeMethod::Volume => {
                        segments
                            .iter_mut()
                            .for_each(|seg| seg.voxelize(space.resolution, bead_radius));
                        segments.sort_by_cached_key(|seg| -> usize {
                            // We can safely unwrap because we just voxelized all segments.
                            seg.voxels().unwrap().count::<true>()
                        });
                        // TODO: Perhaps we can reverse _during_ the sorting operation with some trick?
                        segments.reverse();
                    }
                    RearrangeMethod::Moment => {
                        segments.sort_by_cached_key(|seg| {
                            (seg.structure.moment_of_inertia() * 1e6) as i64
                        });
                        // TODO: Perhaps we can reverse _during_ the sorting operation with some trick?
                        segments.reverse();
                    }

                    RearrangeMethod::BoundingSphere => {
                        segments.sort_by_cached_key(|seg| {
                            (seg.structure.bounding_sphere() * 1e6) as i64
                        });
                        // TODO: Perhaps we can reverse _during_ the sorting operation with some trick?
                        segments.reverse();
                    }
                    // Already taken care of above.
                    RearrangeMethod::None => {
                        unreachable!()
                    }
                }
                eprintln!("Done.");
            }

            segments
        };

        // Output.
        let output = Output {
            title: config.output.title,
            path: args.output,
            topol_includes: config.output.topol_includes,
        };

        let general = General {
            seed,
            max_tries_multiplier,
            max_tries_per_rotation_divisor,
            bead_radius,
        };

        Ok(Self {
            general,
            space,
            segments,
            output,

            rng,
            verbose,
            summary: !args.no_summary,
        })
    }

    pub fn check_rules(&self) -> anyhow::Result<()> {
        let mut checked = Vec::new(); // FIXME: BTreeMap?
        for segment in &self.segments {
            let rules = &segment.rules;
            if rules.is_empty() {
                // No rules to check.
                continue;
            }
            if checked.contains(rules) {
                // Already checked this rule.
                continue;
            }

            // Check the rules.
            let distilled = rules::distill(
                rules,
                self.space.dimensions,
                self.space.resolution,
                &self.space.compartments,
            );
            if !distilled.any::<true>() {
                let name = &segment.name;
                bail!("the rules {rules:?} preclude any placement of segment '{name}'");
            }

            checked.push(rules.clone())
        }

        Ok(())
    }
}

// TODO: This should be part of the config parsing.
fn parse_rule(expr: &RuleExpression) -> Result<Rule, ParseRuleError> {
    match expr {
        RuleExpression::Rule(s) => Rule::from_str(s),
        RuleExpression::Or(exprs) => Ok(Rule::Or(
            exprs.iter().map(parse_rule).collect::<Result<_, _>>()?,
        )),
    }
}
