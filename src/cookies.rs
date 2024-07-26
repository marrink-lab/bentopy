use glam::{BVec3, IVec3, UVec3, Vec3, Vec3A};

use crate::placement::{index_3d, iter_3d};
use crate::structure::{BoxVecsExtension, Structure};
use crate::PeriodicMode;

/// A simple spatial data structure that is used to accelerate solvent-structure collision checks.
///
/// Each _cookie_ represents one solvent box, and is constructed in such a manner that all
/// positions that are in or sufficiently close to the solvent beads represented by a cookie are
/// stored in the cookie.
pub struct Cookies {
    dimensions: UVec3,
    cookie_size: Vec3,
    cookies: Box<[Box<[Vec3]>]>,
}

impl Cookies {
    /// Creates a new [`Cookies`].
    ///
    /// Note that the box size for the `structure` are treated as the dimensions when considering
    /// the treatment of atoms according to the [`PeriodicMode`].
    pub fn new(
        structure: &Structure,
        cookie_size: Vec3,
        dimensions: UVec3,
        cutoff: f32,
        mode: PeriodicMode,
    ) -> Self {
        let n_cells = dimensions.x as usize * dimensions.y as usize * dimensions.z as usize;
        let mut cookies = vec![Vec::new(); n_cells];

        // Place the bead's position in the appropriate cookie.
        // Like raisins or pieces of chocolate. But in our case the cookies aren't round but
        // cuboid. Look, I don't even remember why I picked this name.
        let box_dimensions = structure.boxvecs.as_vec3();
        for (i, bead) in structure.atoms.iter().enumerate() {
            let mut pos = bead.position;
            match mode {
                PeriodicMode::Periodic => pos = pos.rem_euclid(box_dimensions),
                PeriodicMode::Ignore => {
                    if pos.cmplt(Vec3::ZERO).any() || pos.cmpgt(box_dimensions).any() {
                        continue;
                    }
                }
                PeriodicMode::Deny => {
                    if pos.cmplt(Vec3::ZERO).any() || pos.cmpgt(box_dimensions).any() {
                        eprintln!("ERROR: Atom {i} with position {pos} lies outside the box {box_dimensions}.");
                        // FIXME: It's rather gross to me to exit from inside this function. A nice
                        // improvement could be to bubble it up as an Err to the caller.
                        std::process::exit(2);
                    }
                }
            }
            // Note that we can safely cast to a UVec3 here, since beads that fall outside the box
            // dimensions have already been translated inside or skipped.
            assert!(pos.is_negative_bitmask() == 0);
            let cell_pos = (pos / cookie_size).floor().as_uvec3();
            let idx = index_3d(cell_pos, dimensions);
            cookies[idx].push(pos);
        }

        // Now, we have to convolve the contents of the cookies to their neighbors. Otherwise, a
        // neighboring bead may be too close to some solvent bead but go unnoticed in the collision
        // check.
        let idimensions = dimensions.as_ivec3();
        let da = Vec3A::from(box_dimensions);
        // NOTE: If we are ever _really_ pressed for memory, we can win quite a bit here. Cloning
        // the cookies here saves some time gluing the neighbors onto the original positions. But,
        // there are ways of removing this clone.
        let mut convolved_cookies = cookies.clone();
        for cell_pos in iter_3d(dimensions) {
            // This closure returns true if the provided position is in cutoff-range of the present
            // cookie.
            let is_in_range = {
                let start = Vec3A::from(cookie_size * cell_pos.as_vec3() - cutoff);
                let end = Vec3A::from(cookie_size * (cell_pos + 1).as_vec3() + cutoff);
                move |v: &Vec3A| v.cmpge(start).all() && v.cmple(end).all()
            };
            let idx = index_3d(cell_pos, dimensions);

            for offset in NEIGHBORS {
                let neighbor_pos = cell_pos.as_ivec3() + offset;
                // Note that these jumps are mutually exclusive.
                let backward_jumps = neighbor_pos.cmplt(IVec3::ZERO);
                let forward_jumps = neighbor_pos.cmpge(idimensions);
                // Sanity check of their mutual exclusivity.
                assert_eq!(backward_jumps & forward_jumps, BVec3::FALSE);

                if !(backward_jumps | forward_jumps).any() {
                    // Default case, no complicated stuff with periodicity to consider.
                    // TODO: Take only the positions that are closer than `cutoff` to the
                    // central cell.
                    let neighbor_pos = neighbor_pos.as_uvec3();
                    let neighbor_idx = index_3d(neighbor_pos, dimensions);
                    let neighbor_content = cookies[neighbor_idx]
                        .iter()
                        .map(|&v| Vec3A::from(v))
                        .filter(is_in_range)
                        .map(|v| Vec3::from(v));
                    convolved_cookies[idx].extend(neighbor_content);
                } else {
                    fn apply(b: BVec3, v: Vec3A) -> Vec3A {
                        UVec3::new(b.x as u32, b.y as u32, b.z as u32).as_vec3a() * v
                    }

                    let wrapped_neighbor_pos = neighbor_pos.rem_euclid(idimensions).as_uvec3();
                    let wrapped_neighbor_idx = index_3d(wrapped_neighbor_pos, dimensions);
                    assert_ne!(idx, wrapped_neighbor_idx);
                    let wrapped_neighbor_content = cookies[wrapped_neighbor_idx]
                        .iter()
                        .map(|&v| {
                            let translation = apply(backward_jumps, -da) + apply(forward_jumps, da);
                            Vec3A::from(v) + translation
                        })
                        .filter(is_in_range)
                        .map(|v| Vec3::from(v));
                    convolved_cookies[idx].extend(wrapped_neighbor_content);
                }
            }
        }

        Self {
            dimensions,
            cookie_size,
            cookies: convolved_cookies
                .into_iter()
                .map(|cookie| cookie.into_boxed_slice())
                .collect(),
        }
    }

    pub fn get(&self, cell_pos: UVec3) -> Option<&[Vec3]> {
        let idx = index_3d(cell_pos, self.dimensions);
        self.cookies.get(idx).map(|cookie| &cookie[..])
    }

    /// Calculate the offset of a cookie at some position from the system origin.
    pub fn offset(&self, cell_pos: UVec3) -> Vec3 {
        cell_pos.as_vec3() * self.cookie_size()
    }

    pub fn cookie_size(&self) -> Vec3 {
        self.cookie_size
    }
}

const NEIGHBORS: [IVec3; 26] = [
    IVec3::new(-1, -1, -1),
    IVec3::new(0, -1, -1),
    IVec3::new(1, -1, -1),
    IVec3::new(-1, 0, -1),
    IVec3::new(0, 0, -1),
    IVec3::new(1, 0, -1),
    IVec3::new(-1, 1, -1),
    IVec3::new(0, 1, -1),
    IVec3::new(1, 1, -1),
    IVec3::new(-1, -1, 0),
    IVec3::new(0, -1, 0),
    IVec3::new(1, -1, 0),
    IVec3::new(-1, 0, 0),
    IVec3::new(1, 0, 0),
    IVec3::new(-1, 1, 0),
    IVec3::new(0, 1, 0),
    IVec3::new(1, 1, 0),
    IVec3::new(-1, -1, 1),
    IVec3::new(0, -1, 1),
    IVec3::new(1, -1, 1),
    IVec3::new(-1, 0, 1),
    IVec3::new(0, 0, 1),
    IVec3::new(1, 0, 1),
    IVec3::new(-1, 1, 1),
    IVec3::new(0, 1, 1),
    IVec3::new(1, 1, 1),
];
