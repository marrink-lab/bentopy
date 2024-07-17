use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::PathBuf;
use std::str::FromStr;

use eightyseven::reader::ReadGro;
use eightyseven::writer::WriteGro;
use glam::{EulerRot, Mat3, Quat, Vec3};
use rand::SeedableRng;

use crate::args::Args;
use crate::config::{
    Configuration, Mask as ConfigMask, RuleExpression, Shape as ConfigShape, TopolIncludes,
};
use crate::mask::{distance_mask, Dimensions, Mask};
use crate::rules::{self, ParseRuleError, Rule};
use crate::session::{Locations, Session};
use crate::structure::Structure;

const ORDER: EulerRot = EulerRot::XYZ;

pub type CompartmentID = String;
pub type Compartments = Vec<Compartment>;
pub type Size = [f32; 3];
pub type Rotation = Mat3;
pub type Voxels = Mask;
pub type Rng = rand::rngs::StdRng; // TODO: Is this the fastest out there?

impl Mask {
    pub(crate) fn create_from_shape(shape: ConfigShape, dimensions: Dimensions) -> Self {
        let [w, h, d] = dimensions.map(|v| v as usize);
        let r = (dimensions.into_iter().min().unwrap() as usize) / 2;
        let c = r as isize;
        let r2 = r.pow(2);

        // TODO: We can use some nice iterator tricks to avoid allocating a big Vec<bool> here.
        // TODO: Profile to see whether we can help the inlining along.
        let mut cells = Vec::with_capacity(w * h * d);
        for z in 0..d as isize {
            for y in 0..h as isize {
                for x in 0..w as isize {
                    let cell = match shape {
                        ConfigShape::Spherical => {
                            // TODO: Profile and inspect asm to see whether this inlines well.
                            ((x - c).pow(2) as usize)
                                + ((y - c).pow(2) as usize)
                                + ((z - c).pow(2) as usize)
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

    fn load_from_path(path: PathBuf) -> io::Result<Self> {
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
        let mut distance_masks = self.distance_masks.borrow_mut();
        distance_masks
            .entry(key)
            .or_insert_with(|| distance_mask(&self.mask, distance));
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

impl State {
    pub fn new(args: Args, config: Configuration) -> io::Result<Self> {
        let verbose = args.verbose;
        let dimensions = config
            .space
            .size
            .map(|d| (d / config.space.resolution) as u64);
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
                            Mask::create_from_shape(shape, dimensions)
                        }
                        ConfigMask::Voxels { path } => {
                            if verbose {
                                eprintln!("\tLoading mask from {path:?}...");
                            }
                            Mask::load_from_path(path)?
                        }
                    },
                    distance_masks: Default::default(),
                })
            })
            .collect::<io::Result<_>>()?;
        let space = Space {
            size: config.space.size,
            dimensions,
            resolution: config.space.resolution,
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
                    let structure = load_molecule(&seg.path)?;
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
                        name: seg.name,
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
            if args.rearrange {
                eprint!("Rearranging segments... ");
                segments
                    .iter_mut()
                    .for_each(|seg| seg.voxelize(space.resolution, bead_radius));
                segments
                    .sort_by_cached_key(|seg| -> usize { seg.voxels().unwrap().count::<true>() });
                // TODO: Perhaps we can reverse _during_ the sorting operation with some trick?
                segments.reverse();
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
