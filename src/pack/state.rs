use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, bail};
use eightyseven::reader::ReadGro;
use eightyseven::writer::WriteGro;
use glam::{EulerRot, Mat3, Quat, U64Vec3, UVec3, Vec3};
use rand::{Rng as _, RngCore, SeedableRng};

use bentopy::core::config::Axes;
pub use bentopy::core::config::CompartmentID;
use bentopy::core::config::legacy::{
    CombinationExpression, Compartment as ConfigCompartment, Config, Mask as ConfigMask, Quantity,
    RuleExpression, Shape as ConfigShape, TopolIncludes,
};

use crate::args::{Args, RearrangeMethod};
use crate::mask::{Dimensions, Mask, distance_mask_grow};
use crate::placement::{Batch, Placement};
use crate::rules::{self, ParseRuleError, Rule};
use crate::session::{Locations, Session};
use crate::structure::Structure;
use crate::voxelize::voxelize;
use crate::{CLEAR_LINE, Summary};

const ORDER: EulerRot = EulerRot::XYZ;

pub type Compartments = Vec<Compartment>;
pub type Size = [f32; 3];
pub type Rotation = Mat3;
pub type Voxels = Mask;
pub type Rng = rand::rngs::StdRng; // TODO: Is this the fastest out there?

impl Mask {
    /// All dimensions in terms of voxels.
    pub(crate) fn create_from_shape(
        shape: ConfigShape,
        dimensions: Dimensions,
        center: Option<UVec3>,
        radius: Option<u32>,
    ) -> Self {
        let [w, h, d] = dimensions.map(|v| v as usize);
        let min_dim = dimensions.into_iter().min().unwrap() as u32; // Has non-zero length.
        let r = radius.unwrap_or(min_dim / 2);
        assert!(
            min_dim >= r,
            "the provided radius ({r} voxels) must not exceed the space's smallest dimension ({min_dim} voxels)"
        );
        let c = center.unwrap_or(UVec3::splat(r as u32)).as_ivec3();
        assert!(U64Vec3::from_array(dimensions).cmpge(c.as_u64vec3()).all());
        let r2 = r.pow(2);

        // TODO: We can use some nice iterator tricks to avoid allocating a big Vec<bool> here.
        // TODO: Profile to see whether we can help the inlining along.
        let mut cells = Vec::with_capacity(w * h * d);
        for z in 0..d as i32 {
            for y in 0..h as i32 {
                for x in 0..w as i32 {
                    let cell = match shape {
                        ConfigShape::Spherical => {
                            // TODO: Profile and inspect asm to see whether this inlines well.
                            ((x - c.x).pow(2) as u32)
                                + ((y - c.y).pow(2) as u32)
                                + ((z - c.z).pow(2) as u32)
                                >= r2
                        }
                        ConfigShape::Cuboid | ConfigShape::None => false,
                    };
                    cells.push(cell);
                }
            }
        }
        let cells = cells.into_boxed_slice();

        Self::from_cells(dimensions, &cells)
    }

    fn load_from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let mut npz = npyz::npz::NpzArchive::open(path)?;
        // FIXME: This error could be handled more gracefully.
        let first_name = npz
            .array_names()
            .next()
            .context("There should be at least one array in the mask voxel map")?
            .to_string();
        let array = npz.by_name(&first_name)?.unwrap(); // We just asserted the name exists.
        let Ok(size): Result<[u64; 3], _> = array.shape().to_vec().try_into() else {
            let shape = array.shape();
            bail!("a voxel map must have three dimensions, found {shape:?}");
        };
        if size.iter().any(|d| d == &0) {
            bail!("a voxel map must have non-zero sized dimensions, found {size:?}")
        }

        let order = array.order();
        let mut cells: Box<[bool]> = array.into_vec()?.into_boxed_slice();
        cells.iter_mut().for_each(|cell| *cell = !(*cell));
        let cells = match order {
            npyz::Order::C => {
                // We need to re-order the cells to get them into the expected Fortran ordering.
                let mut reordered = Vec::with_capacity(cells.len());
                reordered.resize(cells.len(), false);
                let mut cells = IntoIterator::into_iter(cells);
                let mut reordered = reordered.into_boxed_slice();
                let [w, h, d] = size;
                assert_eq!(cells.len() as u64, w * h * d);
                for x in 0..w {
                    for y in 0..h {
                        for z in 0..d {
                            // We know that cells and reordered have the same size, and we have
                            // asserted that `size` corresponds to the same length.
                            reordered[(x + y * w + z * w * h) as usize] = cells.next().unwrap();
                        }
                    }
                }
                reordered
            }
            npyz::Order::Fortran => cells,
        };

        Ok(Self::from_cells(size, &cells))
    }
}

type DistanceMasksKey = u64;

pub struct Compartment {
    pub id: CompartmentID,
    pub mask: Mask,
    // We make use of some internal mutability, here, through the RefCell. This allows us to lazily
    // set up the distance masks, such that they are created only once we need them.
    distance_masks: RefCell<HashMap<DistanceMasksKey, Mask>>,
}

impl Compartment {
    /// Get or create and return a cloned distance mask for some voxel distance.
    ///
    /// Note that the distance is in terms of voxels, not in nm.
    pub fn get_distance_mask(&self, distance: f32) -> Mask {
        let key = distance as u64;
        {
            // In this scope, we mutably borrow the distance_masks. At the end of the scope, the
            // mutable borrow is dropped, and we are okay to perform a non-mutable borrow after it.
            let mut distance_masks = self.distance_masks.borrow_mut();
            distance_masks
                .entry(key)
                .or_insert_with(|| distance_mask_grow(&self.mask, distance as u64));
        }
        // We can safely get at the `key` because if it was empty, we have inserted it above.
        // We clone here because we can't guarantee that the value is not moved by a resize of
        // the HashMap.
        self.distance_masks.borrow().get(&key).unwrap().clone()
    }
}

pub struct Space {
    pub size: Size,
    pub dimensions: Dimensions,
    pub resolution: f32,
    pub compartments: Compartments,
    pub periodic: bool,

    pub(crate) global_background: Mask,
    pub(crate) session_background: Mask,
    /// The previous session's compartment IDs are used to see if a renewal of locations is
    /// necessary between subsequent segment placements.
    ///
    /// When set to `None`, a renewal of locations is due for the next session, regardless of the
    /// previous session's compartment IDs.
    previous_compartments: Option<HashSet<CompartmentID>>,
    previous_rules: Option<Vec<Rule>>,
}

impl Space {
    pub fn enter_session<'s>(
        &'s mut self,
        compartment_ids: impl IntoIterator<Item = CompartmentID>,
        rules: impl IntoIterator<Item = Rule>,
        locations: &'s mut Locations,
        quantity: Quantity,
    ) -> Session<'s> {
        let compartment_ids = HashSet::from_iter(compartment_ids);
        let rules = Vec::from_iter(rules);

        // TODO: Consider caching this volume like we do for the same compartments below.
        //       The volume can just be associated with a set of previous compartments.
        let target = quantity.bake(|| self.volume(&compartment_ids));

        // Set up a new session background if necessary.
        // Otherwise, leave the session background and locations alone. The session background
        // will stay exactly the same, since it was already set up for this set of
        // compartments. The locations are likely still valid.
        let same_previous_compartments = self
            .previous_compartments
            .as_ref()
            .is_some_and(|prev| prev == &compartment_ids);
        let same_previous_rules = self
            .previous_rules
            .as_ref()
            .is_some_and(|prev| prev == &rules);
        if !same_previous_compartments || !same_previous_rules {
            // Clone the global background, which has all structures stamped onto it.
            self.session_background = self.global_background.clone();

            // Apply the compartments to the background.
            if let Some(merge) = self
                .compartments
                .iter()
                .filter(|comp| compartment_ids.contains(&comp.id))
                .map(|comp| comp.mask.clone())
                .reduce(|mut acc, mask| {
                    acc.merge_mask(&mask);
                    acc
                })
            {
                self.session_background.apply_mask(&merge);
            }
            self.previous_compartments = Some(compartment_ids);

            // Apply the rules to the background.
            let rule_mask =
                rules::distill(&rules, self.dimensions, self.resolution, &self.compartments);
            self.session_background |= !rule_mask;
            self.previous_rules = Some(rules);

            // We must renew the locations as well, based on the newly masked session background.
            locations.renew(self.session_background.linear_indices_where::<false>());
        }

        Session::new(self, locations, target)
    }

    /// Determine the free voxel volume for the specified compartments.
    fn volume(&self, compartment_ids: &HashSet<CompartmentID>) -> f64 {
        let free_voxels = self
            .compartments
            .iter()
            .filter(|comp| compartment_ids.contains(&comp.id))
            .map(|comp| comp.mask.clone())
            .reduce(|mut acc, mask| {
                acc.merge_mask(&mask);
                acc
            })
            .map(|mask| mask.count::<false>())
            .unwrap_or(0);

        let voxel_volume = (self.resolution as f64).powi(3);
        free_voxels as f64 * voxel_volume
    }
}

pub struct Segment {
    pub name: String,
    pub tag: Option<String>,
    pub quantity: Quantity,
    pub compartments: Vec<CompartmentID>,
    pub rules: Vec<Rule>,
    pub path: PathBuf,
    pub rotation_axes: Axes,
    structure: Structure,
    /// The initial rotation of the structure must be applied before the random rotation.
    initial_rotation: Rotation,
    /// Invariant: This rotation must satisfy the constraints set by the `rotation_axes` field by
    /// construction.
    rotation: Rotation,
    voxels: Option<Voxels>,
}

impl Segment {
    /// Set a new rotation the [`Segment`].
    ///
    /// This invalidates the voxelization.
    ///
    /// The internal `rotation_axes` are taken into account when storing the rotation, such that a
    /// rotation stored in a `Segment` is always internally consistent.
    pub fn set_rotation(&mut self, rotation: Rotation) {
        // FIXME: Assert it's a well-formed rotation?
        // TODO: This seems slightly hacky, since we are converting the rotation between different
        // formats a couple of times. It should be fine---we just lose an unimportant bit of
        // accuracy on a random rotation---but perhaps it is more wise to store the rotation
        // internally as a quaternion and only convert it to Mat3 when writing to the placement
        // list.
        let axes = &self.rotation_axes;
        let (ax, ay, az) = Quat::from_mat3(&rotation).to_euler(ORDER);
        let (ax, ay, az) = (
            if axes.x { ax } else { 0.0 },
            if axes.y { ay } else { 0.0 },
            if axes.z { az } else { 0.0 },
        );
        self.rotation = Mat3::from_euler(ORDER, ax, ay, az);
        self.voxels = None;
    }

    /// Get the correctly formed rotation of this [`Segment`].
    pub fn rotation(&self) -> Rotation {
        self.rotation * self.initial_rotation
    }

    /// Voxelize this [`Segment`] according to its current rotation.
    ///
    /// The voxelization can be accessed through [`Segment::voxels`].
    pub fn voxelize(&mut self, resolution: f32, radius: f32) {
        self.voxels = Some(voxelize(
            &self.structure,
            self.rotation(),
            resolution,
            radius,
        ));
    }

    /// If available, return a reference to the voxels that represent this [`Segment`].
    pub fn voxels(&self) -> Option<&Voxels> {
        self.voxels.as_ref()
    }
}

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
    pub fn new(args: Args, config: Config) -> anyhow::Result<Self> {
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
            .partition(ConfigCompartment::is_predefined);
        let mut compartments: Vec<Compartment> = predefined
            .into_iter()
            .map(|comp| -> anyhow::Result<_> {
                let mask = match comp.mask {
                    ConfigMask::Shape(shape) => {
                        if verbose {
                            eprintln!("\tConstructing a {shape} mask...");
                        }
                        Mask::create_from_shape(shape, dimensions, None, None)
                    }
                    ConfigMask::Analytical {
                        shape,
                        center,
                        radius,
                    } => {
                        if verbose {
                            eprintln!("\tConstructing a {shape} mask...");
                        }
                        let center = center.map(|c| (Vec3::from_array(c) / resolution).as_uvec3());
                        let radius = radius.map(|r| (r / resolution) as u32);
                        Mask::create_from_shape(shape, dimensions, center, radius)
                    }
                    ConfigMask::Voxels { path } => {
                        if verbose {
                            eprintln!("\tLoading mask from {path:?}...");
                        }
                        Mask::load_from_path(&path)
                            .with_context(|| format!("Failed to load mask {path:?}"))?
                    }
                    ConfigMask::Combination(_) => {
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
            let ConfigMask::Combination(ces) = combination.mask else {
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

    pub fn pack(&mut self, log: &mut impl io::Write) -> anyhow::Result<(Vec<Placement>, Summary)> {
        let start = std::time::Instant::now();
        let mut locations = Locations::new();
        let mut placements = Vec::new();
        let mut summary = Summary::new();

        let state = self;
        let n_segments = state.segments.len();
        for (i, segment) in state.segments.iter_mut().enumerate() {
            if segment.quantity.is_zero() {
                write!(
                    log,
                    "{prefix}({i:>3}/{n_segments}) Skipping attempt to pack 0 instances of segment '{name}'.{suffix}",
                    prefix = if state.verbose { "" } else { CLEAR_LINE },
                    i = i + 1,
                    name = segment.name,
                    suffix = if state.verbose { "\n" } else { "" }
                )?;
                continue;
            }

            write!(
                log,
                "{prefix}({i:>3}/{n_segments}) Attempting to pack {target:>5} instances of segment '{name}'.{suffix}",
                prefix = if state.verbose { "" } else { CLEAR_LINE },
                i = i + 1,
                target = segment.quantity,
                name = segment.name,
                suffix = if state.verbose { "\n" } else { "" }
            )?;

            // Prepare the session.
            let start_session = std::time::Instant::now();
            let mut session = state.space.enter_session(
                // FIXME: This cloned stuff does not sit right with me.
                segment.compartments.iter().cloned(),
                segment.rules.iter().cloned(),
                &mut locations,
                segment.quantity,
            );

            // Set up the placement record for this segment.
            let mut placement = Placement::new(
                segment.name.clone(),
                segment.tag.clone(),
                segment.path.clone(),
            );

            let mut hits = 0;
            let mut tries = 0; // The number of unsuccessful tries.
            let mut tries_per_rotation = 0;
            // The maximum number of unsuccessful tries.
            let max_tries = state.general.max_tries_multiplier * session.target() as u64;
            // The number of
            let max_tries_per_rotation = max_tries / state.general.max_tries_per_rotation_divisor;
            let mut batch_positions = Vec::new();
            'placement: while hits < session.target() {
                if tries >= max_tries {
                    if state.verbose {
                        writeln!(log, "Exiting after {tries} unsuccessful tries.")?;
                    }
                    break 'placement; // "When you try your best, but you don't succeed."
                }

                // TODO: This can become more efficient for successive placement failures.
                // TODO: Also, this should become a method on Segment.
                let resolution = session.resolution();
                let voxels = match segment.voxels() {
                    Some(voxels) => voxels,
                    None => {
                        segment.voxelize(resolution, state.general.bead_radius);
                        segment.voxels().unwrap() // We just voxelized.
                    }
                };

                // FIXME: Do all this math with glam::U64Vec?
                let [bx, by, bz] = session.dimensions();
                let [sx, sy, sz] = voxels.dimensions();
                if sx > bx || sy > by || sz > bz {
                    tries += 1;

                    // What about another rotation?
                    let rotation = Mat3::from_quat(state.rng.random());
                    segment.set_rotation(rotation);
                    tries_per_rotation = 0; // Reset the counter.
                    if state.verbose {
                        eprintln!("\tNew rotation. The previous rotation would never have fit.")
                    }

                    continue 'placement; // Segment voxels size exceeds the size of the background.
                }
                let [maxx, maxy, maxz] = [bx - sx, by - sy, bz - sz];

                // Pick a random location.
                let position = 'location: loop {
                    let candidate = match session.pop_location(&mut state.rng) {
                        Some(location) => location,
                        None => {
                            // We've gone through all the locations.
                            // TODO: This is rather unlikely, but must be dealt with correctly. I think
                            // this may require a nice custom type that does some internal mutabilty
                            // trickery.
                            // For now, we'll just break here.
                            break 'placement;
                        }
                    };

                    // Convert the linear index location to a spatial index.
                    let position = session.position(candidate).unwrap();

                    // If we are not packing in a periodic manner, we skip candidate positions that
                    // would necessitate periodic behavior. When a location passes this point, it is
                    // valid for non-periodic placement.
                    if !session.periodic() {
                        // Check if the segment will fit in the background, given this position.
                        let [x, y, z] = position;
                        if x >= maxx || y >= maxy || z >= maxz {
                            continue 'location;
                        }
                    }

                    // All seems good, so we continue to see if this position will be accepted.
                    break position;
                };

                // Reject if it would cause a collision.
                if !session.check_collisions(voxels, position) {
                    tries += 1;
                    tries_per_rotation += 1;
                    if tries_per_rotation >= max_tries_per_rotation {
                        let rotation = Mat3::from_quat(state.rng.random());
                        segment.set_rotation(rotation);
                        tries_per_rotation = 0; // Reset the counter.
                    }
                    continue 'placement; // Reject due to collision.
                }

                // We found a good spot. Stomp on those stamps!
                session.stamp(voxels, position);
                // Transform the location to nm.
                batch_positions.push(position.map(|v| v as f32 * session.resolution()));

                // Let's write out the batch and rotate the segment again.
                // TODO: Perhaps we'll need a little transpose here.
                let batch = Batch::new(segment.rotation(), batch_positions.clone());
                placement.push(batch);
                batch_positions.clear();

                let rotation = Mat3::from_quat(state.rng.random());
                segment.set_rotation(rotation);
                tries_per_rotation = 0; // Reset the counter.

                hits += 1;
            }
            let duration = start_session.elapsed().as_secs_f64();

            if state.verbose {
                let total = start.elapsed().as_secs();
                writeln!(
                    log,
                    "                      Packed {hits:>5} instances in {duration:6.3} s. [{total} s] <{:.4} of max_tries, {:.4} of target>",
                    tries as f32 / max_tries as f32,
                    tries as f32 / session.target() as f32
                )?;
            }

            // Save the batches.
            summary.push(segment.name.clone(), session.target(), hits, duration);
            placements.push(placement);
        }

        if !state.verbose {
            writeln!(log)?; // Go down one line to prevent overwriting the last line.
        }

        Ok((placements, summary))
    }
}

fn parse_rule(expr: &RuleExpression) -> Result<Rule, ParseRuleError> {
    match expr {
        RuleExpression::Rule(s) => Rule::from_str(s),
        RuleExpression::Or(exprs) => Ok(Rule::Or(
            exprs.iter().map(parse_rule).collect::<Result<_, _>>()?,
        )),
    }
}

/// Load a [`Structure`] from a structure file.
///
/// This function will return an error if the structure is empty. Downstream functions may assume
/// that a [`Structure`] has at least one atom.
fn load_molecule<P: AsRef<std::path::Path> + std::fmt::Debug>(
    path: P,
) -> anyhow::Result<Structure> {
    let file = std::fs::File::open(&path)?;

    let structure = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gro") => Structure::read_from_file(file)
            .with_context(|| format!("Failed to parse gro file {path:?}"))?,
        Some("pdb") => Structure::read_from_pdb_file(file)
            .with_context(|| format!("Failed to parse PDB file {path:?}"))?,
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

mod combinations {
    use super::*;

    /// Apply some operation over the mask resulting from the first [`CombinationExpression`] in
    /// `exprs` for each subsequent `expr`.
    ///
    /// # Invariants
    ///
    /// Assumes `exprs` has at least one entry.
    fn fold<F: Fn(&mut Mask, Mask)>(
        exprs: &[CombinationExpression],
        compartments: &[Compartment],
        op: F,
    ) -> anyhow::Result<Mask> {
        // Invariant: Assumes `exprs` has at least one entry.
        let (first, exprs) = exprs.split_first().unwrap();
        let mut mask = execute(first, compartments)?;
        // Now apply `op` between the first mask and those resulting from the subsequent `exprs`.
        for expr in exprs {
            let src = execute(expr, compartments)?;
            op(&mut mask, src);
        }

        Ok(mask)
    }

    pub fn execute(
        expr: &CombinationExpression,
        compartments: &[Compartment],
    ) -> anyhow::Result<Mask> {
        let output = match expr {
            // TODO: This is where a Cow may be great!
            CombinationExpression::Id(id) => {
                let compartment = compartments
                    .iter()
                    .find(|c| &c.id == id)
                    .ok_or(anyhow::anyhow!("mask with id {id:?} not (yet) defined"))?;
                compartment.mask.clone()
            }
            CombinationExpression::Not(expr) => !execute(expr, compartments)?,
            CombinationExpression::Union(exprs) => {
                fold(exprs, compartments, |mask, src| *mask &= src)?
            }
            CombinationExpression::Intersect(exprs) => {
                fold(exprs, compartments, |mask, src| *mask |= src)?
            }
        };

        Ok(output)
    }
}
