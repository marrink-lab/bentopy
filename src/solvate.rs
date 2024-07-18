use eightyseven::writer::WriteGro;
use glam::UVec3;
use glam::Vec3;

use crate::placement::Cookies;
use crate::placement::PlaceMap;
use crate::structure::{BoxVecsExtension, Structure};

/// Solvate a [`Structure`] with a template solvent box.
pub fn solvate(
    structure: &mut Structure,
    solvent: &Structure,
    cutoff: f32,
    center: bool,
    // posions: u32,
    // negions: u32,
) {
    // Determine how many copies of the solvent cell are required to fill the input box for each
    // direction. If the solvent box does not fit precisely (viz., the input box size is not an
    // integer multiple of the solvent box size), the size will be overshot. This means that in
    // many cases, the size of the output box may be greater than that of the input structure.
    let d = structure.boxvecs.as_vec3() / solvent.boxvecs.as_vec3();
    let remainder = d.fract() * solvent.boxvecs.as_vec3();
    let dimensions = d.ceil().as_uvec3();
    if center {
        // Add half of the overshot to the positions in the structure to place it in the center of
        // the final box.
        let offset = remainder * 0.5;
        for bead in &mut structure.atoms {
            bead.position += offset;
        }
    }
    let cutoff2 = cutoff.powi(2); // Square to avoid taking the square root of the distances.

    let mut placemap = PlaceMap::new(solvent, dimensions);

    // Cut the input structure into cell-sized cookies that are ordered in the same manner as the
    // solvent placement map is.
    let cookies = Cookies::new(structure, solvent.boxvecs.as_vec3(), dimensions);

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
                // FIXME: Par?
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
                        // let mut bead = *solvent_bead;
                        // bead.pos = solvent_pos;
                        // structure.beads.push(bead)
                    }
                }
            }
        }
    }

    for z_cell in 0..dimensions.z {
        for y_cell in 0..dimensions.y {
            for x_cell in 0..dimensions.x {
                let cell_pos = UVec3::new(x_cell, y_cell, z_cell);
                let translation = cookies.offset(cell_pos);
                let placement = placemap.get(cell_pos).unwrap();
                let solvent_positions = placemap
                    .solvent
                    .atoms() // FIXME: Par?
                    .zip(placement)
                    .filter_map(|(&sb, occupied)| {
                        if occupied {
                            None
                        } else {
                            let mut sb = sb; // Copy the solvent bead.
                            sb.position += translation;
                            Some(sb)
                        }
                    });
                structure.atoms.extend(solvent_positions)
            }
        }
    }

    match &mut structure.boxvecs {
        eightyseven::structure::BoxVecs::Short(three) => {
            *three = (Vec3::from_array(*three) + remainder).to_array()
        }
        eightyseven::structure::BoxVecs::Full(_) => todo!(),
    };
}
