use eightyseven::writer::WriteGro;
use glam::{UVec3, Vec3};

use crate::placement::Cookies;
use crate::placement::PlaceMap;
use crate::structure::{BoxVecsExtension, Structure};
use crate::{BoundaryMode, PeriodicMode};

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
    eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());

    // Deal with the boundary conditions.
    match boundary_mode {
        BoundaryMode::Cut => {
            let scutoff = cutoff;
            let scutoff2 = scutoff.powi(2);

            eprint!("\tCutting solvent to fit box... ");
            let start = std::time::Instant::now();
            for axis in [Axis::X, Axis::Y, Axis::Z] {
                let max = dimensions - 1;
                let box_dimensions = structure.boxvecs.as_vec3();
                let (main, asz, bsz) = match axis {
                    Axis::X => (Vec3::X, dimensions.y, dimensions.z),
                    Axis::Y => (Vec3::Y, dimensions.x, dimensions.z),
                    Axis::Z => (Vec3::Z, dimensions.x, dimensions.y),
                };

                // Set up a 3×3 layer of atoms as a crown on the box dimensions at the `this` box.
                let limit = axis.select(box_dimensions);
                let crown: Box<[Vec3]> = axis
                    .crown(box_dimensions)
                    .map(|translation| {
                        placemap
                            .solvent
                            .atoms()
                            .map(move |a| a.position + translation)
                            // We are only interested in positions that lie at the periodic
                            // interface. Anything beyond the solvent cutoff distance from that
                            // interface is out of reach.
                            .filter(|&position| axis.select(position) <= limit + scutoff)
                    })
                    .into_iter()
                    .flatten()
                    .collect();

                let translation = cookies.offset(main.as_uvec3() * max);
                let this = placemap
                    .solvent
                    .atoms()
                    .map(|a| a.position)
                    .map(|position| position + translation);

                // We'll make a list of solvent atoms we need to reject by index.
                let rejected: Box<[_]> = this
                    .enumerate()
                    .filter_map(|(idx, bead)| {
                        let outside = axis.select(bead) > limit;
                        let too_close = crown
                            .iter()
                            .any(|&crown_bead| bead.distance_squared(crown_bead) < scutoff2);
                        if outside || too_close {
                            Some(idx)
                        } else {
                            None
                        }
                    })
                    .collect();

                // Apply the rejections to all solvent boxes on this face.
                for b in 0..bsz {
                    for a in 0..asz {
                        let pos = match axis {
                            Axis::X => UVec3::new(max.x, a, b),
                            Axis::Y => UVec3::new(a, max.y, b),
                            Axis::Z => UVec3::new(a, b, max.z),
                        };
                        let mut pm = placemap.get_mut(pos).unwrap();
                        for &idx in &rejected {
                            pm.set_occupied(idx)
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

#[derive(Clone, Copy)]
enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    const fn select(&self, v: Vec3) -> f32 {
        v.to_array()[*self as usize]
    }

    const fn roll(&self, mut v: Vec3) -> Vec3 {
        let mut n = 0;
        while n < *self as usize {
            v = Vec3 {
                x: v.z,
                y: v.x,
                z: v.y,
            };
            n += 1;
        }
        v
    }

    const fn crown(&self, box_dimensions: Vec3) -> [Vec3; 9] {
        let d = self.select(box_dimensions);
        let mut crown = [
            Vec3::new(d, -1.0, -1.0),
            Vec3::new(d, 0.0, -1.0),
            Vec3::new(d, 1.0, -1.0),
            Vec3::new(d, -1.0, 0.0),
            Vec3::new(d, 0.0, 0.0),
            Vec3::new(d, 1.0, 0.0),
            Vec3::new(d, -1.0, 1.0),
            Vec3::new(d, 0.0, 1.0),
            Vec3::new(d, 1.0, 1.0),
        ];
        let mut i = 0;
        while i < crown.len() {
            crown[i] = self.roll(crown[i]);
            i += 1;
        }
        crown
    }
}
