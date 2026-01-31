use std::path::Path;

use anyhow::{Context, Result, bail};
use bentopy::core::config::{self, legacy::Shape as LegacyConfigShape};
use glam::{U64Vec3, UVec3, Vec3};

use crate::mask::{Dimensions, Mask};

impl Mask {
    /// All dimensions in terms of voxels.
    pub(crate) fn legacy_create_from_shape(
        shape: LegacyConfigShape,
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
        let c = center.unwrap_or(UVec3::splat(r)).as_ivec3();
        assert!(U64Vec3::from_array(dimensions).cmpge(c.as_u64vec3()).all());
        let r2 = r.pow(2);

        // TODO: We can use some nice iterator tricks to avoid allocating a big Vec<bool> here.
        // TODO: Profile to see whether we can help the inlining along.
        let mut cells = Vec::with_capacity(w * h * d);
        for z in 0..d as i32 {
            for y in 0..h as i32 {
                for x in 0..w as i32 {
                    let cell = match shape {
                        LegacyConfigShape::Spherical => {
                            // TODO: Profile and inspect asm to see whether this inlines well.
                            ((x - c.x).pow(2) as u32)
                                + ((y - c.y).pow(2) as u32)
                                + ((z - c.z).pow(2) as u32)
                                >= r2
                        }
                        LegacyConfigShape::Cuboid | LegacyConfigShape::None => false,
                    };
                    cells.push(cell);
                }
            }
        }
        let cells = cells.into_boxed_slice();

        Self::from_cells(dimensions, &cells)
    }

    // TODO: Still shitty version. Will be better with a rewrite.
    pub(crate) fn create_cuboid(
        dimensions: [u64; 3],
        resolution: f32,
        start: config::Anchor,
        end: config::Anchor,
    ) -> Self {
        let [w, h, d] = dimensions.map(|v| v as usize);
        let dimensions = U64Vec3::from_array(dimensions);

        let start = match start {
            config::Anchor::Start => U64Vec3::ZERO,
            config::Anchor::Center => dimensions / 2,
            config::Anchor::End => dimensions,
            config::Anchor::Point(point) => (Vec3::from(point) / resolution).as_u64vec3(),
        };
        let end = match end {
            config::Anchor::Start => U64Vec3::ZERO,
            config::Anchor::Center => dimensions / 2,
            config::Anchor::End => dimensions,
            config::Anchor::Point(point) => (Vec3::from(point) / resolution).as_u64vec3(),
        };
        // TODO: Consider if we want to do this here, or make the correct ordering an invariant.
        let (start, end) = (U64Vec3::min(start, end), U64Vec3::max(start, end));

        assert!(
            dimensions.cmpge(start).all(),
            "start ({start}) must be within the space dimensions ({dimensions})"
        );
        assert!(
            dimensions.cmpge(end).all(),
            "end ({end}) must be within the space dimensions ({dimensions})"
        );

        // TODO: We can use some nice iterator tricks to avoid allocating a big Vec<bool> here.
        // TODO: Profile to see whether we can help the inlining along.
        let mut cells = vec![true; w * h * d];
        for z in 0..d {
            for y in 0..h {
                for x in 0..w {
                    let idx = x + w * y + w * h * z;
                    let free = (start.x..end.x).contains(&(x as u64))
                        && (start.y..end.y).contains(&(y as u64))
                        && (start.z..end.z).contains(&(z as u64));
                    cells[idx] = !free;
                }
            }
        }
        let cells = cells.into_boxed_slice();

        Self::from_cells(dimensions.into(), &cells)

        // TODO: Go with something along these lines. Likely faster and more memory efficient.
        // let mut mask = Self::new(dimensions);
        // mask.apply_function(|[x,y,z]| );
    }

    pub(crate) fn create_from_shape(
        dimensions: [u64; 3],
        resolution: f32,
        shape: config::Shape,
    ) -> Mask {
        // TODO: The current implementation just giving the helm over to legacy_create_from_shape
        // is very sad and will be corrected.
        match shape {
            config::Shape::Sphere { center, radius } => {
                let center = match center {
                    config::Anchor::Center => (U64Vec3::from(dimensions) / 2).as_uvec3(),
                    config::Anchor::Point(point) => (Vec3::from(point) / resolution).as_uvec3(),
                    config::Anchor::Start => UVec3::ZERO,
                    config::Anchor::End => U64Vec3::from(dimensions).as_uvec3(),
                };
                let radius = (radius / resolution) as u32;
                Self::legacy_create_from_shape(
                    LegacyConfigShape::Spherical,
                    dimensions,
                    Some(center),
                    Some(radius),
                )
            }
            config::Shape::Cuboid { start, end } => {
                Self::create_cuboid(dimensions, resolution, start, end)
            }
        }
        // LegacyConfigMask::Analytical {
        //     shape,
        //     center,
        //     radius,
        // } => {
        //     if verbose {
        //         eprintln!("\tConstructing a {shape} mask...");
        //     }
        //     let center = center.map(|c| (Vec3::from_array(c) / resolution).as_uvec3());
        //     let radius = radius.map(|r| (r / resolution) as u32);
        //     Mask::legacy_create_from_shape(shape, dimensions, center, radius)
        // }
    }

    pub(crate) fn load_from_path<P: AsRef<Path>>(
        path: P,
        expected_dimensions: Dimensions,
    ) -> Result<Self> {
        let mut npz = npyz::npz::NpzArchive::open(path)?;
        // FIXME: This error could be handled more gracefully.
        let first_name = npz
            .array_names()
            .next()
            .context("There should be at least one array in the mask voxel map")?
            .to_string();
        let array = npz.by_name(&first_name)?.unwrap(); // We just asserted the name exists.
        let Ok(dimensions): Result<[u64; 3], _> = array.shape().to_vec().try_into() else {
            let shape = array.shape();
            bail!("a voxel map must have three dimensions, found {shape:?}");
        };
        if dimensions != expected_dimensions {
            bail!(
                "the dimensions of the voxel map do not match those defined in the configuration, \
                    expected {expected_dimensions:?}, found {dimensions:?}"
            )
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
                let [w, h, d] = dimensions;
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

        Ok(Self::from_cells(dimensions, &cells))
    }
}
