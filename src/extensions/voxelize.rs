use glam::{IVec3, U64Vec3, Vec3};
use numpy::{
    ndarray::{Array2, Array3},
    Ix2, Ix3, PyArray, ToPyArray,
};
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyfunction]
/// Return the voxels corresponding to a point cloud.
///
/// Take in a (3, _n_)-shaped `points` array where _n_ >= 1.
///
/// The `resolution` is a value in nm which will become the final voxel side length.
pub(crate) fn py_voxelize<'py>(
    py: Python<'py>,
    points: &PyArray<f64, Ix2>,
    resolution: f64,
    radius: f64,
) -> PyResult<&'py PyArray<f32, Ix3>> {
    voxelize(py, points, resolution, radius)
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

fn voxelize<'py>(
    py: Python<'py>,
    points: &PyArray<f64, Ix2>,
    resolution: f64,
    radius: f64,
) -> PyResult<&'py PyArray<f32, Ix3>> {
    // Verify that the shape of the points matches our expectations.
    let npoints = match points.shape() {
        &[3, npoints] if npoints > 0 => npoints,
        &[3, npoints @ 0] => {
            return Err(PyValueError::new_err(format!(
                "points array must have at least one point, found {npoints}"
            )));
        }
        _ => {
            return Err(PyValueError::new_err(format!(
                "points array should have the shape (3, n), found {:?}",
                points.shape()
            )));
        }
    };

    // First, a very naive method.
    let points = points.readonly();
    let points = points.as_array();
    let min = points
        .rows()
        .into_iter()
        // FIXME: Think about the implications of total_cmp for our case here.
        .map(|d| d.into_iter().copied().min_by(|a, b| a.total_cmp(b)))
        .collect::<Option<Vec<_>>>()
        .expect("there is at least one point");
    let min = Array2::from_shape_vec((3, 1), min).expect("a point has three values") - radius;
    // Translate the points such that they are all above (0, 0, 0).
    let points = points.to_owned() - min;
    let scaled_points = points / resolution; // Transform points according to the resolution.

    let ns = neigbors(radius as f32);
    let mut filled_indices = Vec::with_capacity(npoints * (1 + ns.len()));
    for point in scaled_points.columns() {
        let p = Vec3::from_array(std::array::from_fn(|idx| point[idx] as f32));
        filled_indices.push(p.as_u64vec3()); // Push the point itself.
        filled_indices.extend(ns.map(|d| (p + d).as_u64vec3()));
    }

    let voxels_shape = filled_indices
        .iter()
        .copied()
        .reduce(|max, idx| U64Vec3::max(max, idx))
        .expect("there is at least one point")
        .to_array()
        .map(|v| v as usize + 1);
    let mut voxels = Array3::<f32>::zeros(voxels_shape);
    for idx in filled_indices {
        // NOTE: f32 as usize rounds to 0 and gives us 0 for negative values. We want that.
        voxels[idx.to_array().map(|v| v as usize)] = 1.0;
    }

    Ok(voxels.to_pyarray(py))
}

const CELLS: [IVec3; 27] = [
    // The middle cell.
    IVec3::ZERO,
    // The six sides.
    IVec3::X,
    IVec3::Y,
    IVec3::Z,
    IVec3::NEG_X,
    IVec3::NEG_Y,
    IVec3::NEG_Z,
    // The eight diagonal corners.
    IVec3::new(-1, -1, -1),
    IVec3::new(1, -1, -1),
    IVec3::new(-1, 1, -1),
    IVec3::new(1, 1, -1),
    IVec3::new(-1, -1, 1),
    IVec3::new(1, -1, 1),
    IVec3::new(-1, 1, 1),
    IVec3::new(1, 1, 1),
    // The twelve edge neigbors.
    // xy plane.
    IVec3::new(1, 1, 0),
    IVec3::new(-1, 1, 0),
    IVec3::new(1, -1, 0),
    IVec3::new(-1, -1, 0),
    // xz plane.
    IVec3::new(1, 0, 1),
    IVec3::new(-1, 0, 1),
    IVec3::new(1, 0, -1),
    IVec3::new(-1, 0, -1),
    // yz plane.
    IVec3::new(0, 1, 1),
    IVec3::new(0, -1, 1),
    IVec3::new(0, 1, -1),
    IVec3::new(0, -1, -1),
];

struct Territory<T> {
    size: [usize; 3],
    data: Box<[Vec<T>]>,
}

impl<T> Territory<T> {
    /// Create a new [`Territory<T>`] with size `x`, `y`, `z`.
    fn new(x: usize, y: usize, z: usize) -> Self {
        let n = x * y * z;
        let mut data = Vec::with_capacity(n);
        data.resize_with(n, Vec::new);
        Self {
            size: [x, y, z],
            data: data.into_boxed_slice(),
        }
    }

    fn get(&self, idx: [usize; 3]) -> Option<&Vec<T>> {
        if !self.contains(idx) {
            return None;
        }

        self.data.get(self.index(idx))
    }

    fn get_mut(&mut self, idx: [usize; 3]) -> Option<&mut Vec<T>> {
        if !self.contains(idx) {
            return None;
        }

        self.data.get_mut(self.index(idx))
    }

    /// Returns whether the provided `idx` points to a valid location in this [`Territory<T>`].
    fn contains(&self, idx: [usize; 3]) -> bool {
        self.size.into_iter().zip(idx).all(|(s, i)| s > i)
    }

    /// Return the index into the data for a given `[x, y, z]` index.
    ///
    /// The returned index is valid iff the provided `idx` is valid. This can be checked using
    /// [`Territory<T>::contains`].
    fn index(&self, idx: [usize; 3]) -> usize {
        let [x, y, z] = idx;
        let [width, height, _depth] = self.size;
        x + y * width + z * height * width
    }
}
