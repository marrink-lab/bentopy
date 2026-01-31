use eightyseven::writer::WriteGro;
use glam::{U64Vec3, Vec3A};

use crate::state::{Rotation, Voxels};
use crate::structure::Structure;

/// Returns whether a sphere intersects with an axis-aligned, unit-sized cube.
#[inline(always)]
fn overlap(cube_center: Vec3A, sphere_center: Vec3A, sphere_radius_squared: f32) -> bool {
    let c = 0.5; // Half side-length of the cube.
    let rsc = sphere_center - cube_center;
    let l = Vec3A::max(rsc.abs() - c, Vec3A::ZERO);
    let l_squared = l * l;
    l_squared.element_sum() <= sphere_radius_squared
}

/// Conservatively voxelize a [`Structure`].
///
/// Invariant: A [`Structure`] has at least one atom.
pub fn voxelize(structure: &Structure, rotation: Rotation, resolution: f32, radius: f32) -> Voxels {
    // Extract and rotate the points from the structure.
    // And make them aligned (Vec3A) to encourage vectorization as well.
    let mut points: Box<[Vec3A]> =
        structure.atoms().map(|&position| rotation.mul_vec3a(position.into())).collect();

    // Invariant: A Structure has at least one atom.
    let npoints = points.len();
    assert!(npoints >= 1, "a Structure has at least one atom");

    // Invariant: A Structure has at least one atom.
    let min = points.iter().copied().reduce(|a, b| a.min(b)).unwrap();
    // Subtract one bead radius such that the lower bounds are treated correctly.
    let min = min - radius;
    // Translate the points such that they are all above (0, 0, 0).
    points.iter_mut().for_each(|p| *p -= min);
    // Find the bounds of our voxel grid.
    // Invariant: A Structure has at least one atom.
    let max = points.iter().copied().reduce(|a, b| a.max(b)).unwrap();
    let max = max + radius;

    // Set up our voxels.
    let dimensions = (max / resolution).ceil().as_u64vec3();
    let mut voxels = Voxels::new(dimensions.to_array());

    // Transform points according to the resolution.
    let scaled_points = {
        points.iter_mut().for_each(|p| *p /= resolution);
        points
    };
    // Scale the radius for the resolution along with the points.
    let scaled_radius = radius / resolution;
    let scaled_radius_squared = scaled_radius.powi(2);

    // Go through each point and set the voxels its sphere radius occupies.
    for point in scaled_points {
        let lower_bound = point - scaled_radius;
        let upper_bound = point + scaled_radius;
        debug_assert!(Vec3A::cmpgt(lower_bound + 1e-5, Vec3A::ZERO).all());
        debug_assert!(Vec3A::cmplt(upper_bound - 1e-5, max / resolution).all());
        let lower_bound = lower_bound.floor().as_u64vec3();
        let upper_bound = upper_bound.ceil().as_u64vec3();

        for z in lower_bound.z..upper_bound.z {
            for y in lower_bound.y..upper_bound.y {
                for x in lower_bound.x..upper_bound.x {
                    let position = U64Vec3::new(x, y, z);
                    let cube_center = Vec3A::from(position.as_vec3() + 0.5);
                    if overlap(cube_center, point, scaled_radius_squared) {
                        let idx = position.to_array();
                        voxels.set(idx, true);
                    }
                }
            }
        }
    }

    // Invariant: A Structure has at least one atom.
    debug_assert!(
        voxels.any::<true>(),
        "the produced voxel mask contains no occupied voxels, should be impossible"
    );
    voxels
}
