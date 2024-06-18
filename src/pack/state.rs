use std::io;

use eightyseven::{reader::ReadGro, structure::Structure};
use glam::{Mat3, Vec3};
use rand::SeedableRng;

use crate::args::Args;
use crate::config::{Configuration, Mask as ConfigMask, Output, Shape as ConfigShape};

pub type CompartmentID = String;
pub type Size = [f32; 3];
pub type Dimensions = [u64; 3]; // TODO: Make into usize? Rather awkward in places right now.
pub type Rotation = Mat3;
pub type Voxels = Mask; // TODO: This should become a bitmap.
type Rng = rand::rngs::StdRng; // TODO: Is this the fastest out there?

// TODO: Bitmap optimization.
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

    fn get_mut(&mut self, idx: Dimensions) -> Option<&mut bool> {
        let [w, h, _d] = self.dimensions.map(|v| v as usize);
        let [x, y, z] = idx.map(|v| v as usize);
        self.cells.get_mut(x + y * w + z * w * h)
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
                                < r2
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

    fn load_from_path(path: std::path::PathBuf) -> io::Result<Self> {
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
}

pub struct Segment {
    pub name: String,
    pub number: usize,
    pub compartments: Vec<CompartmentID>,
    structure: Structure,
    rotation: Rotation,
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
        };

        let bead_radius = args.bead_radius;
        let segments = {
            let mut segments: Vec<_> = config
                .segments
                .into_iter()
                .map(|seg| -> io::Result<_> {
                    let file = std::fs::File::open(seg.path)?;
                    let structure = Structure::read_from_file(file)?;
                    Ok(Segment {
                        name: seg.name,
                        number: seg.number,
                        compartments: seg.compartments,
                        structure,
                        rotation: Rotation::IDENTITY,
                        voxels: None,
                    })
                })
                .collect::<io::Result<_>>()?;
            if args.rearrange {
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
        .map(|atom| rotation.mul_vec3(Vec3::from_array(atom.position.to_array())))
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
