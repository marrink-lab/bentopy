use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, bail};
use bentopy::core::config::{Config, defaults};
use bentopy::core::placement::{Meta, Placement, PlacementList};
pub(crate) use glam::{EulerRot, Mat3};
use rand::{RngCore, SeedableRng};

pub use bentopy::core::config;

use crate::args::{Args, RearrangeMethod};
use crate::mask::{Mask, distance_mask_grow};
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
    pub topol_includes: Vec<String>,
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
    /// [configuration](Config).
    pub fn new(args: Args, config: Config) -> anyhow::Result<Self> {
        // Read values from the general section of the config. If a command line argument is given,
        // it overwrites the config value. (And the deprecated env vars have the highest priority.)

        // If no seed is provided, use a random seed.
        let seed =
            args.seed.or(config.general.seed).unwrap_or_else(|| Rng::from_os_rng().next_u64());
        let rng = Rng::seed_from_u64(seed);

        let bead_radius =
            args.bead_radius.or(config.general.bead_radius).unwrap_or(defaults::BEAD_RADIUS);

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
            args.max_tries_mult
                .or(config.general.max_tries_mult)
                .unwrap_or(defaults::MAX_TRIES_MULT)
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
                .or(config.general.max_tries_rot_div)
                .unwrap_or(defaults::MAX_TRIES_ROT_DIV)
        };

        let verbose = args.verbose;

        // Space.
        // TODO: Consider if the resolution should be a default, again?
        let resolution =
            config.space.resolution.context("No resolution was specified in the input file")?
                as f32;
        let size =
            config.space.dimensions.context("No dimensions were specified in the input file")?;
        // The dimensions from the config is the real-space size of the box. Here, we treat the
        // word dimensions as being the size in terms of voxels. Bit annoying.
        // TODO: Reconsider this wording.
        let dimensions = size.map(|d| (d / resolution) as u64);

        eprintln!("Setting up compartments...");
        let (predefined, combinations): (Vec<_>, Vec<_>) =
            config.compartments.into_iter().partition(config::Compartment::is_predefined);
        let mut compartments: Vec<Compartment> = predefined
            .into_iter()
            .map(|comp| -> anyhow::Result<_> {
                let mask = match comp.mask {
                    config::Mask::Voxels(path) => {
                        if verbose {
                            eprintln!("\tLoading mask from {path:?}...");
                        }
                        Mask::load_from_path(&path, dimensions)
                            .with_context(|| format!("Failed to load mask {path:?}"))?
                    }
                    config::Mask::All => {
                        if verbose {
                            eprintln!("\tConstructing a full space mask...");
                        }
                        let shape = config::Shape::Cuboid {
                            start: config::Anchor::Start,
                            end: config::Anchor::End,
                        };
                        Mask::create_from_shape(dimensions, resolution, shape)
                    }
                    config::Mask::Shape(shape) => {
                        if verbose {
                            eprintln!("\tConstructing a {shape} mask...");
                        }
                        Mask::create_from_shape(dimensions, resolution, shape)
                    }
                    config::Mask::Limits(expr) => {
                        if verbose {
                            eprintln!("\tConstructing a mask from limits...");
                        }
                        compartment::distill_limits(&expr, dimensions, resolution as f64)
                    }
                    // We partitioned the list, so these variants are not present.
                    config::Mask::Within { .. }
                    | config::Mask::Around { .. }
                    | config::Mask::Combination(_) => unreachable!(),
                };

                Ok(Compartment { id: comp.id, mask })
            })
            .collect::<anyhow::Result<_>>()?;
        for combination in combinations {
            if verbose {
                eprintln!("\tApplying a compartment combination...");
            }

            let baked = match combination.mask {
                // Within some distance of the mask for compartment `id`. This is inclusive.
                config::Mask::Within { distance, id } => {
                    let compartment = compartments
                        .iter()
                        .find(|c| c.id == id)
                        .ok_or(anyhow::anyhow!("Mask with id {id:?} not (yet) defined"))?;
                    let mask = &compartment.mask;
                    let voxel_distance = (distance / resolution) as u64;
                    Compartment {
                        id: combination.id,
                        mask: distance_mask_grow(mask, voxel_distance),
                    }
                }
                // Mask within some distance of compartment `id` excluding that source mask.
                // This is exclusive. Around is Within followed by a not and on the source mask.
                config::Mask::Around { distance, id } => {
                    let compartment = compartments
                        .iter()
                        .find(|c| c.id == id)
                        .ok_or(anyhow::anyhow!("Mask with id {id:?} not (yet) defined"))?;
                    let source = &compartment.mask;
                    let voxel_distance = (distance / resolution) as u64;
                    // First, the same step as for within.
                    let mut mask = distance_mask_grow(source, voxel_distance);
                    // Now we need to take away the source mask from the within mask.
                    // FIXME: I'm still so frustrated with the inversed mask decision. Ugh past me.
                    // That's why the logic here looks so cursed. Practically, OR means AND for the
                    // way I designed the mask datatype.
                    mask |= !source.clone();
                    Compartment { id: combination.id, mask }
                }
                config::Mask::Combination(expr) => Compartment {
                    id: combination.id,
                    mask: combinations::execute(&expr, &compartments)?,
                },

                // We partitioned the list, so these variants are not present.
                config::Mask::All
                | config::Mask::Voxels(_)
                | config::Mask::Shape(_)
                | config::Mask::Limits(_) => unreachable!(),
            };

            compartments.push(baked);
        }
        let space = Space {
            size,
            dimensions,
            resolution,
            compartments,
            periodic: config.space.periodic.unwrap_or(defaults::PERIODIC),

            global_background: Mask::new(dimensions),
            session_background: Mask::new(dimensions),
            previous_compartments: None,
        };

        // Segments.
        eprintln!("Loading segment structures...");
        let segments = {
            let constraints: HashMap<&String, &config::Rule> =
                config.constraints.iter().map(|c| (&c.id, &c.rule)).collect();
            let mut segments: Vec<_> = config
                .segments
                .into_iter()
                .map(|seg| -> Result<_, _> {
                    let path = seg.path;
                    if verbose {
                        eprintln!("\tLoading {path:?}...");
                    }
                    let name = seg.name;
                    let tag = seg.tag;
                    match tag.as_ref().map(String::len) {
                        Some(0) => eprintln!("WARNING: The tag for segment '{name}' is empty."),
                        Some(6.. ) => eprintln!("WARNING: The tag for segment '{name}' is longer than 5 characters, and may be truncated when the placement list is rendered."),
                        _ => {} // Nothing to warn about.
                    }
                    // TODO: Use the segment.name() method here? Or rather, it's better named
                    //       future version ;)
                    let structure = load_molecule(&path).with_context(|| format!("Failed to open the structure file for segment '{name}' at {path:?}"))?;

                    // If there are any, there must be only one. We don't enforce that here though.
                    let rotation_axes = if let Some(id) = seg.rules.first() {
                         match constraints
                            .get(id)
                            .ok_or(anyhow::anyhow!("Constraint with id {id:?} is not defined"))? {
                                config::Rule::RotationAxes(axes) => *axes,
                            }
                    } else {
                            config::Axes::default()
                        }
                    ;

                    Ok(Segment {
                        name,
                        tag,
                        quantity: seg.quantity,
                        compartments: seg.compartment_ids,
                        path,
                        rotation_axes,
                        structure,

                        // TODO: Entirely remove initial_rotation from the Segment struct?
                        initial_rotation: Rotation::IDENTITY,
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
                            .for_each(|seg| seg.voxelize(space.resolution, bead_radius as f32));
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
            title: config.general.title.unwrap_or(defaults::TITLE.to_string()),
            path: args.output,
            // TODO: One more example of using a PathBuf here being a poor choice.
            // TODO: Quite some cleanup here, eventually.
            topol_includes: config
                .includes
                .into_iter()
                .map(|p| p.to_string_lossy().to_string())
                .collect(),
        };

        let general = General {
            seed,
            max_tries_multiplier,
            max_tries_per_rotation_divisor,
            bead_radius: bead_radius as f32,
        };

        Ok(Self { general, space, segments, output, rng, verbose, summary: !args.no_summary })
    }

    pub fn check_masks(&self) -> anyhow::Result<()> {
        if self.verbose {
            eprintln!("Checking compartments...")
        }

        let mut checked = Vec::new(); // FIXME: BTreeMap?
        for segment in &self.segments {
            let ids = &segment.compartments;
            let ids_formatted = ids.join(", ");
            let name = &segment.name;
            if ids.is_empty() {
                // No compartments to check?!
                unreachable!("a segment has at least one compartment")
            }
            if checked.contains(ids) {
                // Already checked this rule.
                continue;
            }

            // Check the masks.
            let masks =
                self.space.compartments.iter().filter(|c| ids.contains(&c.id)).map(|c| &c.mask);
            // TODO: Check for correctness.
            let distilled = masks
                .cloned()
                .reduce(|mut acc, m| {
                    acc.merge_mask(&m);
                    acc
                })
                // We know there is at least one compartment.
                // TODO: Don't we already check this upstream?
                .ok_or(anyhow::anyhow!(
                    "No valid compartments ({ids_formatted}) declared for segment '{name}'"
                ))?;

            if self.verbose {
                let n = distilled.count::<false>();
                eprintln!("\t{n:>12} open voxels in compartments '{ids_formatted}'.")
            }

            if !distilled.any::<false>() {
                // Provide some extra debug info about this failing segment's compartments if there
                // is more than one. This is helpful in debugging problems if this error is hit.
                if segment.compartments.len() > 1 {
                    eprintln!("\tIndividual compartments for segment '{name}':");
                    let compartments =
                        self.space.compartments.iter().filter(|c| ids.contains(&c.id));
                    for compartment in compartments {
                        let id = &compartment.id;
                        let n = compartment.mask.count::<false>();
                        eprintln!("\t{n:>12} open voxels in compartment '{id}'");
                    }
                }
                bail!(
                    "The compartments '{ids_formatted}' together preclude any placement of segment '{name}'"
                );
            }

            checked.push(ids.clone())
        }

        if self.verbose {
            let n = checked.len();
            eprintln!("\tAll okay. Checked {n} segment compartment combinations.")
        }

        Ok(())
    }

    pub fn placement_list(&self, placements: impl IntoIterator<Item = Placement>) -> PlacementList {
        let meta = Meta {
            seed: self.general.seed,
            max_tries_mult: self.general.max_tries_multiplier,
            max_tries_per_rotation_divisor: self.general.max_tries_per_rotation_divisor,
            bead_radius: self.general.bead_radius,
        };

        PlacementList {
            title: self.output.title.to_string(),
            size: self.space.size,
            meta: Some(meta),
            topol_includes: self.output.topol_includes.clone(),
            placements: placements.into_iter().collect(),
        }
    }
}
