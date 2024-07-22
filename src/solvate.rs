use eightyseven::writer::WriteGro;
use glam::BVec3;
use glam::UVec3;
use glam::Vec3;

use crate::placement::Cookies;
use crate::placement::PlaceMap;
use crate::structure::{BoxVecsExtension, Structure};
use crate::BoundaryMode;
use crate::PeriodicMode;

/// Solvate a [`Structure`] with a template solvent box.
pub fn solvate<'sol>(
    structure: &mut Structure,
    solvent: &'sol Structure,
    cutoff: f32,
    center: bool,
    boundary_mode: BoundaryMode,
    periodic_mode: PeriodicMode,
) -> PlaceMap<'sol> {
    // Determine how many copies of the solvent cell are required to fill the input box for each
    // direction.
    let d = structure.boxvecs.as_vec3() / solvent.boxvecs.as_vec3();
    let dimensions = d.ceil().as_uvec3();

    // Depending on the desired behavior, modify the structure.
    match boundary_mode {
        BoundaryMode::Cut => {}
        BoundaryMode::Grow => {
            eprintln!("\tGrowing output structure box to fit solvent precisely.");
            // If the solvent box does not fit precisely (viz., the input box size is not an
            // integer multiple of the solvent box size), the size will be overshot. This means that in
            // many cases, the size of the output box may be greater than that of the input structure.
            let remainder = d.fract() * solvent.boxvecs.as_vec3();
            if center {
                // Add half of the overshot to the positions in the structure to place it in the center of
                // the final box.
                let offset = remainder * 0.5;
                for bead in &mut structure.atoms {
                    bead.position += offset;
                }

                eprintln!("\tCentered structure in expanded box.");
            }

            match &mut structure.boxvecs {
                // We don't have to cut off the solvent structures at the boundaries, but we have to
                // change the boundaries to fit.
                eightyseven::structure::BoxVecs::Short(three) => {
                    *three = (Vec3::from_array(*three) + remainder).to_array()
                }
                eightyseven::structure::BoxVecs::Full(_) => todo!(),
            }
        }
    }

    // Cut the input structure into cell-sized cookies that are ordered in the same manner as the
    // solvent placement map is.
    eprint!("\tCutting cookies... ");
    let start = std::time::Instant::now();
    let cookies = Cookies::new(
        structure,
        solvent.boxvecs.as_vec3(),
        dimensions,
        cutoff,
        periodic_mode,
    );
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    // For each location of a solvent box, go through each of the solvent bead positions and see
    // whether it collides with a structure atom. If it does, we jot that down so we don't write
    // that solvent bead out later.
    eprint!("\tChecking solvent collisions... ");
    let start = std::time::Instant::now();
    let mut placemap = PlaceMap::new(solvent, dimensions);
    let cutoff2 = cutoff.powi(2); // Square to avoid taking the square root of the distances.
    for z_cell in 0..dimensions.z {
        for y_cell in 0..dimensions.y {
            for x_cell in 0..dimensions.x {
                let cell_pos = UVec3::new(x_cell, y_cell, z_cell);
                let cookie = cookies.get(cell_pos).unwrap(); // We are sure there is a cookie here.
                if cookie.is_empty() {
                    // If the cookie does not contain any structure beads, we can just skip
                    // measuring any distances.
                    continue;
                }
                let translation = cookies.offset(cell_pos);
                for (idx, solvent_bead) in placemap.solvent.atoms().enumerate() {
                    let solvent_pos = solvent_bead.position + translation;
                    let mut collision = false;
                    for &cookie_bead in cookie {
                        // TODO: Consider applying this translation as a map before the iter over
                        // the solvent beads?
                        // Translate the unit cell solvent bead position to the reference frame of
                        // the bead from the cookie.
                        if solvent_pos.distance_squared(cookie_bead) < cutoff2 {
                            collision = true;
                            break;
                        }
                    }
                    if collision {
                        // This unwrap is safe, since we are sure there is a cell here.
                        let placement = &mut placemap.get_mut(cell_pos).unwrap();
                        placement.set_occupied(idx)
                    }
                }
            }
        }
    }
    eprintln!("Took {:.6} s.", start.elapsed().as_secs_f32());

    // Deal with the boundary conditions.
    match boundary_mode {
        BoundaryMode::Cut => {
            // We will cut the boundaries as follows.
            // For the three outward planes (x, y, z), and for the three outward edges (xy, xz, yz),
            // we can determine the adjustment once and copy it over the entire side or edge.
            // The outer corner (xyz) we take care of with just making the adjustment for that
            // single box.

            eprint!("\tCutting solvent to fit box... ");
            let start = std::time::Instant::now();
            for axes in OUTERS {
                // FIXME: Turn this into a proper type or an extension on MutPlacement or something.
                let mut rejected = Vec::new();

                let uaxes = UVec3::new(axes.x as u32, axes.y as u32, axes.z as u32);
                let vaxes = uaxes.as_vec3();

                let boxvecs = structure.boxvecs.as_vec3();
                let max = dimensions - 1;
                let cell_pos = max * uaxes;
                let translation = cookies.offset(cell_pos);
                // The periodic neigbours.
                let periodic_translation = translation + boxvecs * vaxes;
                let other_side: Box<[_]> = placemap
                    .solvent
                    .atoms()
                    .map(|atom| atom.position)
                    // We filter the beads to make sure that we're not needlessly checking against
                    // beads that are not at the edge of the periodic boundaries.
                    .filter(|&p| (p * vaxes).to_array().into_iter().all(|v| v <= cutoff))
                    .map(|p| p + periodic_translation)
                    .collect();
                for (idx, solvent_bead) in placemap.solvent.atoms().enumerate() {
                    let solvent_pos = solvent_bead.position + translation;

                    // Reject the bead if it lies outside the box.
                    let outside = (solvent_pos.cmpgt(boxvecs) & axes).any();

                    // Check for collisions with the periodic image.
                    let mut collision = false;
                    for &periodic_bead in &other_side {
                        if solvent_pos.distance_squared(periodic_bead) < cutoff2 {
                            collision = true;
                            break;
                        }
                    }

                    // FIXME: Check whether this is compiled as expected, that is, when outside is
                    // true, it is an early out into this branch. In other words, collision is lazy.
                    if outside || collision {
                        rejected.push(idx)
                    }
                }

                // Set the rejected beads to occupied in the placemap.
                // FIXME: I don't like this. I don't like this. I don't like this. See above about
                // `rejected`.
                let [xrange, yrange, zrange] = ranges(dimensions, axes);
                for z in zrange {
                    for y in yrange.clone() {
                        for x in xrange.clone() {
                            let cell_pos = UVec3::new(x, y, z);
                            let mut pm = placemap.get_mut(cell_pos).unwrap();
                            for &idx in &rejected {
                                pm.set_occupied(idx)
                            }
                        }
                    }
                }
            }
            eprintln!("Took {:.6} s.", start.elapsed().as_secs_f32());
        }
        // We already set the box size at the top. Nothing to do here.
        BoundaryMode::Grow => {}
    }

    placemap
}

const OUTERS: [BVec3; 7] = [
    // The planes (x, y, z).
    BVec3::new(true, false, false),
    BVec3::new(false, true, false),
    BVec3::new(false, false, true),
    // The edges (xy, xz, yz).
    BVec3::new(true, true, false),
    BVec3::new(true, false, true),
    BVec3::new(false, true, true),
    // The outer corner (xyz).
    BVec3::new(true, true, true),
];

fn ranges(dimensions: UVec3, axes: BVec3) -> [std::ops::Range<u32>; 3] {
    let m = (dimensions - 1).to_array(); // The maximum index.
    let d = dimensions.to_array(); // A shorter alias.
    [0, 1, 2].map(|i| if axes.test(i) { m[i]..d[i] } else { 0..m[i] })
}
