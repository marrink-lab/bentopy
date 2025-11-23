use std::path::Path;

use anyhow::{Context, Result, bail};
use bentopy::core::config::Shape;
use bentopy::core::config::legacy::Shape as LegacyConfigShape;
use glam::{U64Vec3, UVec3};

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

    /// All dimensions in terms of voxels.
    pub(crate) fn _create_from_shape(_dimensions: Dimensions, _shape: Shape) -> Self {
        todo!()
    }

    pub(crate) fn load_from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
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
