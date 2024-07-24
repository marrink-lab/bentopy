use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use eightyseven::reader::ReadGro;
use eightyseven::writer::WriteGro;
use glam::{EulerRot, Mat3, Quat, U64Vec3, UVec3, Vec3};
use rand::{Rng as _, SeedableRng};

use crate::args::{Args, RearrangeMethod};
use crate::config::{
    Configuration, Mask as ConfigMask, RuleExpression, Shape as ConfigShape, TopolIncludes,
};
use crate::mask::{distance_mask, distance_mask_grow, Dimensions, Mask};
use crate::placement::{Batch, Placement};
use crate::rules::{self, ParseRuleError, Rule};
use crate::session::{Locations, Session};
use crate::structure::Structure;
use crate::{Summary, CLEAR_LINE};

const ORDER: EulerRot = EulerRot::XYZ;

pub type CompartmentID = String;
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
        let min_dim = dimensions.into_iter().min().unwrap() as u32;
        let r = radius.unwrap_or(min_dim / 2);
        assert!(min_dim >= r);
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

    fn load_from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut npz = npyz::npz::NpzArchive::open(path)?;
        // FIXME: This error could be handled more gracefully.
        let first_name = npz
            .array_names()
            .next()
            .expect("there should be at least one array in the mask voxel map")
            .to_string();
        let array = npz.by_name(&first_name)?.unwrap(); // We just asserted the name exists.
        assert_eq!(
            array.shape().len(),
            3,
            "a voxel map must have three dimensions"
        );
        let size: [u64; 3] = array.shape().to_vec().try_into().unwrap(); // We just asserted there are three items.

        let order = array.order();
        let mut cells: Box<[bool]> = array.into_vec()?.into_boxed_slice();
        cells.iter_mut().for_each(|cell| *cell = !(*cell));
        let cells = match order {
            npyz::Order::C => {
                // We need to re-order the cells to get them into the expected Fortran ordering.
                let mut reordered = Vec::with_capacity(cells.len());
                reordered.resize(cells.len(), false);
                let mut cells = cells.iter().copied();
                let mut reordered = reordered.into_boxed_slice();
                let [w, h, d] = size;
                for x in 0..w {
                    for y in 0..h {
                        for z in 0..d {
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

impl From<ParseRuleError> for io::Error {
    fn from(value: ParseRuleError) -> Self {
        io::Error::other(value)
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
    /// Get or create and return a reference to a distance mask for some voxel distance.
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
        target: usize,
    ) -> Session {
        let compartment_ids = HashSet::from_iter(compartment_ids);
        let rules = Vec::from_iter(rules);

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
            for compartment in self
                .compartments
                .iter()
                .filter(|comp| compartment_ids.contains(&comp.id))
            {
                self.session_background.apply_mask(&compartment.mask)
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
}

pub struct Axes {
    x: bool,
    y: bool,
    z: bool,
}

impl Default for Axes {
    fn default() -> Self {
        Self {
            x: true,
            y: true,
            z: true,
        }
    }
}

impl FromStr for Axes {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let len = s.len();
        if len > 3 {
            return Err(format!(
                "an axes string may consist of at most 3 characters, '{s}' has {len} characters"
            ));
        }

        Ok(Self {
            x: s.contains('x'),
            y: s.contains('y'),
            z: s.contains('z'),
        })
    }
}

pub struct Segment {
    pub name: String,
    pub tag: Option<String>,
    pub target: usize,
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

pub struct State {
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,

    pub rng: Rng,
    pub bead_radius: f32,
    pub verbose: bool,
    pub summary: bool,
}

/// Helper function for reporting clear errors when opening something at a path fails.
fn report_opening<P: std::fmt::Debug>(err: impl std::error::Error, path: P) -> io::Error {
    io::Error::other(format!("problem while opening {path:?}: {err}"))
}

impl State {
    pub fn new(args: Args, config: Configuration) -> io::Result<Self> {
        let verbose = args.verbose;
        let dimensions = config
            .space
            .size
            .map(|d| (d / config.space.resolution) as u64);
        let resolution = config.space.resolution;
        eprintln!("Setting up compartments...");
        let compartments = config
            .space
            .compartments
            .into_iter()
            .map(|comp| -> io::Result<_> {
                Ok(Compartment {
                    id: comp.id,
                    mask: match comp.mask {
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
                            let center =
                                center.map(|c| (Vec3::from_array(c) / resolution).as_uvec3());
                            let radius = radius.map(|r| (r / resolution) as u32);
                            Mask::create_from_shape(shape, dimensions, center, radius)
                        }
                        ConfigMask::Voxels { path } => {
                            if verbose {
                                eprintln!("\tLoading mask from {path:?}...");
                            }
                            Mask::load_from_path(&path).map_err(|err| report_opening(err, path))?
                        }
                    },
                    distance_masks: Default::default(),
                })
            })
            .collect::<io::Result<_>>()?;
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

        let bead_radius = args.bead_radius;
        eprintln!("Loading segment structures...");
        let segments = {
            let mut segments: Vec<_> = config
                .segments
                .into_iter()
                .map(|seg| -> io::Result<_> {
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
                    let structure = load_molecule(&seg.path).map_err(|err| report_opening(err, &seg.path))?;
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
                        target: seg.number,
                        compartments: seg.compartments,
                        rules,
                        path: seg.path,
                        rotation_axes: seg.rotation_axes,
                        structure,
                        initial_rotation,
                        rotation: Rotation::IDENTITY,
                        voxels: None,
                    })
                })
                .collect::<io::Result<_>>()?;
            if let Some(method) = args.rearrange {
                eprint!("Rearranging segments according to the {method:?} method... ");
                match method {
                    RearrangeMethod::Volume => {
                        segments
                            .iter_mut()
                            .for_each(|seg| seg.voxelize(space.resolution, bead_radius));
                        segments.sort_by_cached_key(|seg| -> usize {
                            seg.voxels().unwrap().count::<true>()
                        });
                        // TODO: Perhaps we can reverse _during_ the sorting operation with some trick?
                        segments.reverse();
                    }
                    RearrangeMethod::MomentOfInertia => {
                        segments.sort_by_cached_key(|seg| {
                            (seg.structure.moment_of_inertia() * 1e6) as i64
                        });
                        // TODO: Perhaps we can reverse _during_ the sorting operation with some trick?
                        segments.reverse();
                    }
                }
                eprintln!("Done.");
            }
            segments
        };

        let output = Output {
            title: config.output.title,
            path: args.output,
            topol_includes: config.output.topol_includes,
        };

        let rng = if let Some(seed) = args.seed {
            Rng::seed_from_u64(seed)
        } else {
            Rng::from_entropy()
        };

        Ok(Self {
            space,
            segments,
            output,

            rng,
            bead_radius,
            verbose,
            summary: !args.no_summary,
        })
    }

    pub fn check_rules(&self) -> Result<(), String> {
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
                return Err(format!(
                    "the rules {rules:?} preclude any placement of segment '{name}'"
                ));
            }

            checked.push(rules.clone())
        }

        Ok(())
    }

    pub fn pack(
        &mut self,
        log: &mut impl io::Write,
    ) -> Result<(Vec<Placement>, Summary), io::Error> {
        let start = std::time::Instant::now();
        let mut locations = Locations::new();
        let mut placements = Vec::new();
        let mut summary = Summary::new();

        let max_tries_multiplier = std::env::var("BENTOPY_TRIES")
            .map(|s| {
                s.parse()
                    .expect("max tries multiplier should be a valid unsigned integer")
            })
            .inspect(|n| eprintln!("\tMax tries multiplier set to {n}."))
            .unwrap_or(1000);
        let max_tries_per_rotation_divisor = std::env::var("BENTOPY_ROT_DIV")
            .map(|s| {
                s.parse()
                    .expect("divisor should be a valid unsigned integer")
            })
            .inspect(|n| eprintln!("\tMax tries per rotation divisor set to {n}."))
            .unwrap_or(100);

        let state = self;
        let n_segments = state.segments.len();
        for (i, segment) in state.segments.iter_mut().enumerate() {
            if segment.target == 0 {
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
                target = segment.target,
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
                segment.target,
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
            let max_tries = max_tries_multiplier * segment.target; // The number of unsuccessful tries.
            let max_tries_per_rotation = max_tries / max_tries_per_rotation_divisor;
            let mut batch_positions = Vec::new();
            'placement: while hits < segment.target {
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
                        segment.voxelize(resolution, state.bead_radius);
                        segment.voxels().unwrap()
                    }
                };

                // FIXME: Do all this math with glam::U64Vec?
                let [bx, by, bz] = session.dimensions();
                let [sx, sy, sz] = voxels.dimensions();
                if sx > bx || sy > by || sz > bz {
                    tries += 1;

                    // What about another rotation?
                    let rotation = Mat3::from_quat(state.rng.gen());
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
                        let rotation = Mat3::from_quat(state.rng.gen());
                        segment.set_rotation(rotation);

                        if state.verbose {
                            eprintln!("\tNew rotation. Exceeded the maximum number of tries per rotation ({tries_per_rotation} @ {tries} tries).")
                        }

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

                let rotation = Mat3::from_quat(state.rng.gen());
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
                    tries as f32 / segment.target as f32
                )?;
            }

            // Save the batches.
            summary.push(segment.name.clone(), segment.target, hits, duration);
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
fn load_molecule<P: AsRef<std::path::Path> + std::fmt::Debug>(path: P) -> io::Result<Structure> {
    let file = std::fs::File::open(&path)?;

    let structure = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gro") => Structure::read_from_file(file)?,
        Some("pdb") => Structure::read_from_pdb_file(file)?,
        None | Some(_) => {
            eprintln!("WARNING: Assuming {path:?} is a pdb file.");
            Structure::read_from_pdb_file(file)?
        }
    };

    Ok(structure)
}

/// Calculate the translations to check whether a bead would partially occupy a neigboring voxel
/// for a given `radius`.
///
/// This does not include the identity translation.
fn neigbors(radius: f32) -> [Vec3; 26] {
    // Calculate the distance for which the magnitude of the vector is equal to `radius`.
    let t2d = f32::sqrt((1.0 / 2.0) * radius.powi(2));
    let t3d = f32::sqrt((1.0 / 3.0) * radius.powi(2));

    [
        // The six sides.
        Vec3::X * radius,
        Vec3::Y * radius,
        Vec3::Z * radius,
        Vec3::NEG_X * radius,
        Vec3::NEG_Y * radius,
        Vec3::NEG_Z * radius,
        // The eight diagonal corners.
        Vec3::new(-t3d, -t3d, -t3d),
        Vec3::new(t3d, -t3d, -t3d),
        Vec3::new(-t3d, t3d, -t3d),
        Vec3::new(t3d, t3d, -t3d),
        Vec3::new(-t3d, -t3d, t3d),
        Vec3::new(t3d, -t3d, t3d),
        Vec3::new(-t3d, t3d, t3d),
        Vec3::new(t3d, t3d, t3d),
        // The twelve edge neigbors.
        // xy plane.
        Vec3::new(t2d, t2d, 0.0),
        Vec3::new(-t2d, t2d, 0.0),
        Vec3::new(t2d, -t2d, 0.0),
        Vec3::new(-t2d, -t2d, 0.0),
        // xz plane.
        Vec3::new(t2d, 0.0, t2d),
        Vec3::new(-t2d, 0.0, t2d),
        Vec3::new(t2d, 0.0, -t2d),
        Vec3::new(-t2d, 0.0, -t2d),
        // yz plane.
        Vec3::new(0.0, t2d, t2d),
        Vec3::new(0.0, -t2d, t2d),
        Vec3::new(0.0, t2d, -t2d),
        Vec3::new(0.0, -t2d, -t2d),
    ]
}

fn voxelize(structure: &Structure, rotation: Rotation, resolution: f32, radius: f32) -> Voxels {
    // TODO: Consider whether this rotations system is correct. I tend to flip things around.
    // Extract and rotate the points from the structure.
    let mut points: Box<[_]> = structure
        .atoms()
        .map(|&position| rotation.mul_vec3(position))
        .collect();
    let npoints = points.len();
    // FIXME: Think about the implications of total_cmp for our case here.
    let min = points
        .iter()
        .copied()
        .reduce(|a, b| a.min(b))
        .expect("there is at least one point");
    // Subtract one bead radius such that the lower bounds are treated correctly.
    let min = min - radius;
    // Translate the points such that they are all above (0, 0, 0).
    points.iter_mut().for_each(|p| *p -= min);

    // Transform points according to the resolution.
    let scaled_points = points.iter().map(|p| *p / resolution);
    // Scale the radius for the resolution along with the points.
    let scaled_radius = radius / resolution;

    let ns = neigbors(scaled_radius);
    let mut filled_indices = Vec::with_capacity(npoints * (1 + ns.len()));
    for point in scaled_points {
        // TODO: This as_u64vec3() suddenly seems suspect to me. Investigate whether this actually
        // makes sense and does not produce overly eager voxelizations at the box edges.
        filled_indices.push(point.as_u64vec3()); // Push the point itself.
        filled_indices.extend(ns.map(|d| (point + d).as_u64vec3())); // And the neighbors.
    }

    let voxels_shape = filled_indices
        .iter()
        .copied()
        .reduce(|a, b| a.max(b))
        .expect("there is at least one point")
        .to_array()
        .map(|v| v + 1);
    let mut voxels = Voxels::new(voxels_shape);
    for idx in filled_indices {
        voxels.set(idx.to_array(), true);
    }

    voxels
}
