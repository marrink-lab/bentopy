use eightyseven::writer::WriteGro;
use glam::Vec3;

use crate::state::{Rotation, Voxels};
use crate::structure::Structure;

/// Voxelize a [`Structure`].
///
/// Invariant: A [`Structure`] has at least one atom.
pub fn voxelize(structure: &Structure, rotation: Rotation, resolution: f32, radius: f32) -> Voxels {
    // TODO: Consider whether this rotations system is correct. I tend to flip things around.
    // Extract and rotate the points from the structure.
    let mut points: Box<[_]> = structure
        .atoms()
        .map(|&position| rotation.mul_vec3(position))
        .collect();
    // Invariant: A Structure has at least one atom.
    let npoints = points.len();
    assert!(npoints >= 1, "a Structure has at least one atom");
    // FIXME: Think about the implications of total_cmp for our case here.
    let min = points
        .iter()
        .copied()
        .reduce(|a, b| a.min(b))
        // Invariant: A Structure has at least one atom.
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
        // Invariant: A Structure has at least one atom.
        .expect("there is at least one point")
        .to_array()
        .map(|v| v + 1);
    let mut voxels = Voxels::new(voxels_shape);
    for idx in filled_indices {
        voxels.set(idx.to_array(), true);
    }

    voxels
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
