use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;

use arraystring::typenum::U5;
use arraystring::ArrayString;
use eightyseven::reader::{ParseList, ReadGro};
use glam::{Mat3, U64Vec3, Vec3};
use rand::SeedableRng;

use crate::args::Args;
use crate::config::{Configuration, Mask as ConfigMask, Output, Shape as ConfigShape};

pub type CompartmentID = String;
pub type Size = [f32; 3];
pub type Dimensions = [u64; 3]; // TODO: Make into usize? Rather awkward in places right now.
pub type Position = Dimensions;
pub type Rotation = Mat3;
pub type Voxels = Mask; // TODO: This should become a bitmap.
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

// TODO: Bitmap optimization.
#[derive(Clone)]
pub struct Mask {
    dimensions: Dimensions,
    cells: Box<[bool]>,
}

impl Mask {
    fn new(dimensions: Dimensions) -> Self {
        let [w, h, d] = dimensions.map(|v| v as usize);
        let cells = vec![false; w * h * d].into_boxed_slice();
        Self { dimensions, cells }
    }

    pub const fn dimensions(&self) -> Dimensions {
        self.dimensions
    }

    pub const fn dimensions_usize(&self) -> [usize; 3] {
        let [x, y, z] = self.dimensions;
        [x as usize, y as usize, z as usize]
    }

    // TODO: Impl get and get_mut as Deref<Slice>? Would that be nicer?
    #[inline(always)]
    fn get(&self, idx: Position) -> Option<&bool> {
        let [w, h, _d] = self.dimensions_usize();
        let [x, y, z] = idx.map(|v| v as usize);
        self.cells.get(x + y * w + z * w * h)
    }

    #[inline(always)]
    fn get_mut(&mut self, idx: Position) -> Option<&mut bool> {
        let [w, h, _d] = self.dimensions_usize();
        let [x, y, z] = idx.map(|v| v as usize);
        self.cells.get_mut(x + y * w + z * w * h)
    }

    /// Return a three-dimensional index into this [`Mask`] from a linear index.
    ///
    /// This function makes no guarantees about whether the returned index is actually within the
    /// `Mask`.
    #[inline(always)]
    fn idx(&self, mut i: usize) -> Position {
        let [w, h, _d] = self.dimensions_usize();

        let x = i % w;
        i /= w;
        let y = i % h;
        let z = i / h;

        [x as u64, y as u64, z as u64]
    }

    /// Return an [`Iterator`] over all indices where the the cell is equal to `VALUE`.
    pub fn indices_where<const VALUE: bool>(&self) -> impl Iterator<Item = Position> + '_ {
        self.linear_indices_where::<VALUE>().map(|i| self.idx(i))
    }

    /// Return an [`Iterator`] over all linear indices into the internal `cells` where the the cell
    /// is equal to `VALUE`.
    pub fn linear_indices_where<const VALUE: bool>(&self) -> impl Iterator<Item = usize> + '_ {
        self.cells
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v == VALUE { Some(i) } else { None })
    }

    /// Apply some `mask` onto this [`Mask`].
    ///
    /// This is essentially an `or-assign` operation between the `self` and the provided `mask`.
    fn apply_mask(&mut self, mask: &Mask) {
        assert_eq!(self.dimensions, mask.dimensions);
        assert_eq!(self.cells.len(), mask.cells.len()); // For good measure, so the compiler gets it.
        self.cells
            .iter_mut()
            .zip(mask.cells.iter())
            .for_each(|(s, &m)| *s |= m);
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

        Self { dimensions, cells }
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
        let cells = array.into_vec()?.into_boxed_slice();
        Ok(Self {
            dimensions: size,
            cells,
        })
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
    session_background: Option<Mask>,
}

impl Space {
    pub fn stamp(&mut self, voxels: &Mask, location: Position) {
        let location = U64Vec3::from_array(location);
        for idx in voxels
            .indices_where::<true>()
            .map(|idx| U64Vec3::from_array(idx) + location)
        {
            *self
                .global_background
                .get_mut(idx.to_array())
                .expect("cannot stamp outside of bounds") = true;
            *self
                .session_background
                .as_mut()
                .expect("when stamping, a session background must exist")
                .get_mut(idx.to_array())
                .expect("cannot stamp outside of bounds") = true;
        }
    }

    // TODO: enter_session and exit_session are rather good candidates for a typestate pattern.
    pub fn enter_session(&mut self, compartment_ids: &[CompartmentID]) {
        assert!(
            self.session_background.is_none(),
            "hey, let's not enter a session when one is still ongoing"
        );

        // Clone the global background, which has all structures stamped onto it.
        self.session_background
            .replace(self.global_background.clone());
        let session_background = self.session_background.as_mut().unwrap();

        for compartment in self
            .compartments
            .iter()
            .filter(|comp| compartment_ids.contains(&comp.id))
        {
            session_background.apply_mask(&compartment.mask)
        }
    }

    pub fn exit_session(&mut self) {
        assert!(
            self.session_background.is_some(),
            "hey, let's not exit a session when there is not one going on"
        );

        // Forget the session background.
        self.session_background = None;
    }

    pub fn get_free_locations(&self) -> Vec<Position> {
        let b = self
            .session_background
            .as_ref()
            .expect("the session background must exist");

        let mut lin_idx = 0;
        let mut locations = Vec::new();
        let [w, h, d] = b.dimensions();
        for z in 0..d {
            for y in 0..h {
                for x in 0..w {
                    // TODO: Profile to see if this unchecked get can benefit us here.
                    if !b.cells[lin_idx] {
                        locations.push([x, y, z]);
                    }
                    lin_idx += 1;
                }
            }
        }

        locations
    }

    /// Returns `true` if no collisions are encountered between the [`Space`] and the provided
    /// [`Voxels`] at some `location`.
    ///
    /// If a collision is found, the function exits early with a `false` value.
    pub fn check_collisions(&self, voxels: &Voxels, location: Position) -> bool {
        let occupied_indices = voxels.indices_where::<true>().map(U64Vec3::from_array);
        let location = U64Vec3::from_array(location);
        let occupied_indices = occupied_indices.map(|p| p + location);
        let session_background = self
            .session_background
            .as_ref()
            .expect("the session background must exist");
        for idx in occupied_indices.map(|p| p.to_array()) {
            match session_background.get(idx) {
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

    /// If available, take the voxels that represent this [`Segment`].
    pub fn take_voxels(&mut self) -> Option<Voxels> {
        self.voxels.take()
    }

    // TODO: This is obsolete unless we move on to an internal mutability system.
    /// Return a reference to the voxels that represent this [`Segment`].
    ///
    /// If a voxelization already exists, it is returned. Otherwise, a new voxelization is
    /// determined and returned.
    ///
    /// Note that this function operates under the assumption that a potential previous
    /// voxelization was produced with the same resolution and radius. The responsibility for
    /// upholding this assumption is with the caller.
    pub fn get_voxels(&mut self, resolution: f32, radius: f32) -> &Voxels {
        if self.voxels().is_none() {
            self.voxelize(resolution, radius);
        }

        self.voxels().unwrap()
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
            session_background: None,
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
                segments.sort_by_cached_key(|seg| -> usize {
                    seg.voxels()
                        .unwrap()
                        .cells
                        .iter()
                        .map(|&v| v as usize)
                        .sum()
                });
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
        // NOTE: f32 as usize rounds to 0 and gives us 0 for negative values. We want that.
        // FIXME: Ponder whether we can safely unwrap here.
        *voxels.get_mut(idx.to_array()).unwrap() = true;
    }

    voxels
}
