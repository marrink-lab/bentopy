use glam::{I64Vec3, U64Vec3, Vec3};
use numpy::ndarray::s;
use numpy::{Ix3, PyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::PyList;
use rand::seq::SliceRandom;
use rand::SeedableRng;

#[pyfunction]
/// Place instances of a segment onto a background.
pub(crate) fn py_place<'py>(
    py: Python<'py>,
    global_background: &PyArray<f32, Ix3>,
    session_background: &PyArray<f32, Ix3>,
    segment: &PyArray<f32, Ix3>,
    background_offset: &PyArray1<i64>,
    resolution: f32,
    target: usize,
    max_tries: usize,
    seed: u64,
) -> PyResult<&'py PyList> {
    place(
        py,
        global_background,
        session_background,
        segment,
        background_offset,
        resolution,
        target,
        max_tries,
        seed,
    )
}

fn place<'py>(
    py: Python<'py>,
    global_background: &PyArray<f32, Ix3>,
    session_background: &PyArray<f32, Ix3>,
    segment: &PyArray<f32, Ix3>,
    background_offset: &PyArray1<i64>,
    resolution: f32,
    target: usize,
    max_tries: usize,
    seed: u64,
) -> PyResult<&'py PyList> {
    let mut global_background = global_background.readwrite();
    let mut global_background = global_background.as_array_mut();
    let mut session_background = session_background.readwrite();
    let mut session_background = session_background.as_array_mut();
    let segment_voxels = segment.readonly();
    let segment_voxels = segment_voxels.as_array();

    // Get the indices of the open voxels.
    // Make sure we don't consider cells where a placement onto it would place the segment outside
    // of the background.
    let session_background_shape = {
        let s = session_background.shape();
        U64Vec3::new(s[0] as u64, s[1] as u64, s[2] as u64)
    };
    let segment_voxels_shape = {
        let s = segment_voxels.shape();
        U64Vec3::new(s[0] as u64, s[1] as u64, s[2] as u64)
    };
    let l = U64Vec3::min(session_background_shape, segment_voxels_shape);
    let [lx, ly, lz] = l.to_array().map(|v| v as isize);
    let bg = session_background.slice(s![..-lx, ..-ly, ..-lz]);
    let valid = bg
        .indexed_iter()
        .filter_map(|(idx, &v)| if v == 0.0 { Some(idx) } else { None });
    let mut locations: Vec<_> = valid.collect();
    if locations.is_empty() {
        // If there's no valid spots left, return with an empty list.
        return Ok(PyList::empty(py));
    }

    // Set up a shuffled pile of locations.
    // We do a partial shuffle because the number of locations that will be sampled tends to be
    // very low compared to the total number of available spaces.
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    // TODO: This is a very conservative value. Initial testing showed it can probably be much lower!
    let initial_shuffle_estimate = target + max_tries;
    let (mut shuffled, mut unshuffled) =
        locations.partial_shuffle(&mut rng, initial_shuffle_estimate);
    let mut cursor = 0;

    let mut placements = Vec::new();
    let mut hits = 0;
    let mut tries = 0;
    'placement: while hits < target {
        if tries >= max_tries {
            // Give up.
            break;
        }

        // Take a random location from our shuffled pile.
        let location = match shuffled.get(cursor) {
            Some(s) => s,
            // Tried all open options in the current shuffled batch.
            // If there are unshuffled positions left, shuffle those and move on.
            None if !unshuffled.is_empty() => {
                (shuffled, unshuffled) =
                    unshuffled.partial_shuffle(&mut rng, initial_shuffle_estimate);
                cursor = 0;
                &shuffled[cursor]
            }
            // We're done here.
            None => break,
        };
        cursor += 1;

        // Make sure that this placement does not overlap with occupied voxels.
        let prospect: Vec<_> = segment_voxels
            .indexed_iter()
            .filter_map(|(idx, &v)| if v == 1.0 { Some(idx) } else { None })
            .map(|(x, y, z)| (x + location.0, y + location.1, z + location.2))
            .collect();
        for &index in &prospect {
            match session_background.get(index) {
                Some(0.0) => {} // Fine, the cell is free.
                Some(_) => {
                    // Found a non-zero value, so this location is already occupied.
                    tries += 1;
                    continue 'placement;
                }
                // FIXME: Let's report this for now, even though it really should be impossible.
                None => eprintln!("weird: we should always be in range here! {index:?}"),
            }
        }

        // Stamp the background.
        for &index in &prospect {
            global_background[index] = 1.0;
            session_background[index] = 1.0;
            // FIXME: For checking the bounds here, see above.
        }

        // Store the placement location.
        let location = Vec3::new(location.0 as f32, location.1 as f32, location.2 as f32);
        let offset =
            I64Vec3::from_slice(background_offset.to_owned_array().as_slice().unwrap()).as_vec3();
        let corrected_location = (location + offset) * resolution;
        placements.push(corrected_location.to_array());
        hits += 1;
    }

    eprintln!("tries: {tries},\tmax_tries: {max_tries},\thits: {hits},\ttarget: {target}");

    Ok(PyList::new(py, placements))
}
