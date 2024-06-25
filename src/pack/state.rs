use std::collections::HashSet;
use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;

use arraystring::typenum::U5;
use arraystring::ArrayString;
use eightyseven::reader::{ParseList, ReadGro};
use eightyseven::writer::format_atom_line;
use glam::{Mat3, U64Vec3, Vec3};
use rand::seq::SliceRandom;
use rand::SeedableRng;

use crate::args::Args;
use crate::config::{Configuration, Mask as ConfigMask, Output, Shape as ConfigShape};
use crate::mask::{Dimensions, Mask, Position};
use crate::Location;

pub type CompartmentID = String;
pub type Size = [f32; 3];
pub type Rotation = Mat3;
pub type Voxels = Mask;
type Rng = rand::rngs::StdRng; // TODO: Is this the fastest out there?
type Atom = Vec3;

pub struct Structure {
    atoms: Vec<Atom>,
}

impl ReadGro<Atom> for Structure {
    const PARSE_LIST: ParseList = ParseList {
        resnum: false,
        resname: false,
        atomname: false,
        atomnum: false,
        position: true,
        velocity: false,
    };

    // TODO: Fix this ArrayString bs in eightyseven.
    fn build_atom(
        _resnum: Option<u32>,
        _resname: Option<ArrayString<U5>>,
        _atomname: Option<ArrayString<U5>>,
        _atomnum: Option<u32>,
        position: Option<[f32; 3]>,
        _velocity: Option<[f32; 3]>,
    ) -> Atom {
        Atom::from_array(position.unwrap())
    }

    fn build_structure(
        _title: String,
        atoms: Vec<Atom>,
        _boxvecs: eightyseven::structure::BoxVecs,
    ) -> Self {
        Self { atoms }
    }
}

impl Structure {
    fn read_from_pdb_file(
        file: std::fs::File,
    ) -> Result<Structure, eightyseven::reader::ParseGroError> {
        let reader = BufReader::new(file);

        let mut atoms = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let x = line[30..38].trim().parse().unwrap();
                let y = line[38..46].trim().parse().unwrap();
                let z = line[46..54].trim().parse().unwrap();
                let atom = Atom::new(x, y, z) / 10.0; // Convert from Å to nm.
                atoms.push(atom);
            }
        }

        Ok(Self { atoms })
    }
}

use eightyseven::writer::WriteGro;
impl<'atoms> WriteGro<'atoms, Atom> for Structure {
    fn title(&self) -> String {
        "debug".to_string()
    }

    fn natoms(&self) -> usize {
        self.atoms.len()
    }

    fn atoms(&'atoms self) -> impl Iterator<Item = &'atoms Atom> {
        self.atoms.iter()
    }

    fn boxvecs(&self) -> String {
        "400.0 400.0 400.0".to_string()
    }

    fn format_atom_line(atom: &Atom) -> String {
        format_atom_line(1, "DUMMY", "DUMMY", 2, atom.to_array(), None)
    }
}

impl Mask {
    fn create_from_shape(shape: ConfigShape, dimensions: Dimensions) -> Self {
        let [w, h, d] = dimensions.map(|v| v as usize);
        let r = (dimensions.into_iter().min().unwrap() as usize) / 2;
        let c = r as isize;
        let r2 = r.pow(2);

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
        let size = array.shape().to_vec().try_into().unwrap(); // We just asserted there are three items.
        let mut cells = array.into_vec()?.into_boxed_slice();
        cells
            .iter_mut()
            .for_each(|cell: &mut bool| *cell = !(*cell));
        Ok(Self::from_cells(size, &cells))
    }
}

pub struct Compartment {
    id: CompartmentID,
    mask: Mask,
}

pub struct Space {
    pub size: Size,
    pub dimensions: Dimensions,
    pub resolution: f32,
    pub compartments: Vec<Compartment>,

    global_background: Mask,
    session_background: Mask,
    previous_compartments: Option<HashSet<CompartmentID>>,
}

impl Space {
    // TODO: enter_session and exit_session are rather good candidates for a typestate pattern.
    pub fn enter_session<'s>(
        &'s mut self,
        compartment_ids: impl IntoIterator<Item = CompartmentID>,
        locations: &'s mut Locations,
        target: usize,
    ) -> Session {
        let compartment_ids = HashSet::from_iter(compartment_ids);

        // Set up a new session background if necessary.
        // Otherwise, leave the session background and locations alone. The session background
        // will stay exactly the same, since it was already set up for this set of
        // compartments. The locations are likely still valid.
        if !self
            .previous_compartments
            .as_ref()
            .is_some_and(|prev| prev == &compartment_ids)
        {
            // Clone the global background, which has all structures stamped onto it.
            self.session_background = self.global_background.clone();

            for compartment in self
                .compartments
                .iter()
                .filter(|comp| compartment_ids.contains(&comp.id))
            {
                self.session_background.apply_mask(&compartment.mask)
            }
            self.previous_compartments = Some(compartment_ids);

            // We must renew the locations as well, based on the newly masked session background.
            locations.renew(self.session_background.linear_indices_where::<false>());
        }

        Session {
            inner: self,
            locations,
            target,
        }
    }
}

pub struct Session<'s> {
    inner: &'s mut Space,
    locations: &'s mut Locations,
    target: usize,
}

impl Session<'_> {
    pub fn resolution(&self) -> f32 {
        self.inner.resolution
    }

    pub fn dimensions(&self) -> Dimensions {
        self.inner.dimensions
    }

    /// Returns a location if available.
    ///
    /// This function also takes a reference to a random number generator, since the internal
    /// [`Locations`] type may shuffle its indices on demand.
    ///
    /// If the background for this [`Session`] does not have any free locations left, the session
    /// must be ended. This case is indicated with a `None` return value.
    pub fn pop_location(&mut self, rng: &mut Rng) -> Option<Location> {
        // TODO: Consider having the shuffle guess that is passed to pop be dynamic with the number
        // of placements that have already occurred. We could update that number as we go.
        if let Some(location) = self.locations.pop(rng, self.target) {
            return Some(location);
        }

        // We're out of locations. We'll try and renew once.
        self.locations.renew(
            self.inner
                .session_background
                .linear_indices_where::<false>(),
        );

        // If the result is still `None`, there is nothing we can do anymore.
        self.locations.pop(rng, self.target)
    }

    /// Returns `true` if no collisions are encountered between the [`Space`] and the provided
    /// [`Voxels`] at some `position`.
    ///
    /// If a collision is found, the function exits early with a `false` value.
    pub fn check_collisions(&self, voxels: &Voxels, position: Position) -> bool {
        let occupied_indices = voxels.indices_where::<true>().map(U64Vec3::from_array);
        let position = U64Vec3::from_array(position);

        for idx in occupied_indices.map(|p| (p + position).to_array()) {
            match self.inner.session_background.get(idx) {
                Some(false) => continue,    // Free. Moving on to the next one.
                Some(true) => return false, // Already occupied!
                None => {
                    // Out of bounds!?
                    #[cfg(debug_assertions)]
                    eprintln!("out of bounds this is weird");
                    return false;
                }
            }
        }

        return true;
    }

    pub fn stamp(&mut self, voxels: &Mask, location: Position) {
        let location = U64Vec3::from_array(location);
        for idx in voxels
            .indices_where::<true>()
            .map(|idx| U64Vec3::from_array(idx) + location)
        {
            self.inner.global_background.set(idx.to_array(), true);
            self.inner.session_background.set(idx.to_array(), true);
        }
    }

    pub fn exit_session(&mut self) {
        // Not really doing anything here anymore.
    }

    pub const fn position(&self, location: Location) -> Option<Position> {
        self.inner.session_background.spatial_idx(location)
    }
}

impl<'s> Drop for Session<'s> {
    fn drop(&mut self) {
        self.exit_session();
    }
}

pub struct Locations {
    /// Linear indices into the background map.
    indices: Vec<Location>,
    used: usize,
    threshold: usize,

    /// The number of shuffled elements that may be popped from the end of `indices`.
    cursor: usize,
}

impl Locations {
    /// Threshold of spare capacity at which the `locations` [`Vec`] ought to be shrunk.
    const SHRINK_THRESHOLD: usize = (100 * 1024_usize.pow(2)) / std::mem::size_of::<Location>();
    /// Value that determines the threshold at which the locations will be renewed.
    const RENEW_THRESHOLD: f64 = 0.1;

    pub fn new() -> Self {
        Self {
            indices: Vec::new(),
            used: 0,
            threshold: 0,
            cursor: 0,
        }
    }

    /// Renew the locations from an iterator of locations.
    pub fn renew(&mut self, locations: impl Iterator<Item = Location>) {
        let indices = &mut self.indices;
        indices.clear();
        indices.extend(locations);
        // Shrink the `locations` if the spare capacity is getting out of hand.
        if indices.spare_capacity_mut().len() > Self::SHRINK_THRESHOLD {
            indices.shrink_to_fit();
        }
        self.used = 0;
        // We want to make sure that the threshold is rounded down, such that it eventually could
        // reach 0, indicating the enclosing Session is exhausted. See `pop`.
        self.threshold = (indices.len() as f64 * Self::RENEW_THRESHOLD).floor() as usize
    }

    fn shuffle(&mut self, rng: &mut Rng, shuffle_guess: usize) {
        let (shuffled, _unshuffled) = self.indices.partial_shuffle(rng, shuffle_guess);
        self.cursor = shuffled.len();
    }

    /// Pop a location.
    ///
    /// This function also takes a reference to a random number generator, since this type may
    /// shuffle its indices on demand. The `shuffle_guess` suggests the number of items that ought
    /// to be shuffled in that scenario.
    ///
    /// A `None` value indicates that this [`Locations`] object must be [renewed](`Self::renew`),
    /// not necessarily that all locations have been popped.
    pub fn pop(&mut self, rng: &mut Rng, shuffle_guess: usize) -> Option<Location> {
        if self.used >= self.threshold {
            return None; // Indicate a renewal is needed.
        }
        if self.cursor == 0 {
            self.shuffle(rng, shuffle_guess);
        }
        let location = self.indices.pop()?;
        self.used += 1;
        assert!(self.cursor > 0);
        self.cursor -= 1;
        Some(location)
    }
}

pub struct Segment {
    pub name: String,
    pub target: usize,
    pub compartments: Vec<CompartmentID>,
    pub path: PathBuf,
    structure: Structure,
    pub(crate) rotation: Rotation,
    voxels: Option<Voxels>,
}

impl Segment {
    /// Set a new rotation the [`Segment`].
    ///
    /// This invalidates the voxelization.
    pub fn set_rotation(&mut self, rotation: Rotation) {
        // FIXME: Assert it's a well-formed rotation?
        self.rotation = rotation;
        self.voxels = None;
    }

    /// Voxelize this [`Segment`] according to its current rotation.
    ///
    /// The voxelization can be accessed through [`Segment::voxels`].
    pub fn voxelize(&mut self, resolution: f32, radius: f32) {
        self.voxels = Some(voxelize(&self.structure, self.rotation, resolution, radius));
    }

    /// If available, return a reference to the voxels that represent this [`Segment`].
    pub fn voxels(&self) -> Option<&Voxels> {
        self.voxels.as_ref()
    }
}

pub struct State {
    pub space: Space,
    pub segments: Vec<Segment>,
    pub output: Output,

    pub rng: Rng,
    pub rotations: usize,
    pub bead_radius: f32,
    pub verbose: bool, // TODO: Add verbose flag functionality.
}

impl State {
    pub fn new(args: Args, config: Configuration) -> io::Result<Self> {
        let dimensions = config
            .space
            .size
            .map(|d| (d / config.space.resolution) as u64);
        let space = Space {
            size: config.space.size,
            dimensions,
            resolution: config.space.resolution,
            compartments: config
                .space
                .compartments
                .into_iter()
                .map(|comp| -> io::Result<_> {
                    Ok(Compartment {
                        id: comp.id,
                        mask: match comp.mask {
                            ConfigMask::Shape(shape) => Mask::create_from_shape(shape, dimensions),
                            ConfigMask::Voxels { path } => Mask::load_from_path(path)?,
                        },
                    })
                })
                .collect::<io::Result<_>>()?,

            global_background: Mask::new(dimensions),
            session_background: Mask::new(dimensions),
            previous_compartments: None,
        };

        let bead_radius = args.bead_radius;
        let segments = {
            let mut segments: Vec<_> = config
                .segments
                .into_iter()
                .map(|seg| -> io::Result<_> {
                    let structure = load_molecule(&seg.path)?;
                    Ok(Segment {
                        name: seg.name,
                        target: seg.number,
                        compartments: seg.compartments,
                        path: seg.path,
                        structure,
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

        let rng = if let Some(seed) = args.seed {
            Rng::seed_from_u64(seed)
        } else {
            Rng::from_entropy()
        };

        let output_dir = &config.output.dir;
        if !output_dir.exists() {
            // TODO: Reconsider whether we'll want to eprint from here. Seems a bit crunchy.
            eprintln!("Output directory {output_dir:?} does not exist, yet, and will be created.");
            std::fs::create_dir_all(&output_dir)?;
        }

        Ok(Self {
            space,
            segments,
            output: config.output,

            rng,
            rotations: args.rotations,
            bead_radius,
            verbose: args.verbose, // TODO: Add verbose flag functionality.
        })
    }
}

/// Load a [`Structure`] from a structure file.
fn load_molecule<P: AsRef<std::path::Path> + std::fmt::Debug>(path: P) -> io::Result<Structure> {
    eprintln!("\tLoading {path:?}...");
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
        .atoms
        .iter()
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
