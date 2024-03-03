use glam::{UVec3, Vec3};

use crate::structure::Structure;

/// Return an index into a linear array that represents items on a 3-dimensional grid in a z-major
/// ordering.
fn index_3d(pos: UVec3, dimensions: UVec3) -> usize {
    // Casts to usize here serve to avoid wrapping when indexing _very_ large sizes. You never know!
    pos.x as usize
        + pos.y as usize * dimensions.x as usize
        + pos.z as usize * dimensions.x as usize * dimensions.y as usize
}

pub struct PlaceMap<'s> {
    pub solvent: &'s Structure,
    dimensions: UVec3,
    // Invariant: length = (solvent.beads.len() / 8).ceil() * dimensions.product()
    // Value is 0 if the place of this solvent bead is not otherwise occupied.
    // Conversely, if the value of a bit is 1, a solvent bead cannot be placed there.
    placements: Box<[u8]>,
}

impl<'s> PlaceMap<'s> {
    pub fn new(solvent: &'s Structure, size: UVec3) -> Self {
        let n_cells = size.x as usize * size.y as usize * size.z as usize;
        let n_bytes = solvent.beads.len().div_ceil(8);
        let placements = vec![0x00; n_bytes * n_cells].into_boxed_slice();
        Self {
            dimensions: size,
            solvent,
            placements,
        }
    }

    pub fn get_mut(&mut self, pos: UVec3) -> Option<PlacementMut> {
        if pos.x >= self.dimensions.x || pos.y >= self.dimensions.y || pos.z >= self.dimensions.z {
            return None;
        }

        // TODO: Is it wasteful to 'compute' this here? Or does it reduce down well?
        let n_beads = self.solvent.beads.len();
        let n_bytes = n_beads.div_ceil(8);
        let idx = index_3d(pos, self.dimensions);
        // This unwrap should be safe, since we checked at the start of the function.
        let bytes = self.placements.chunks_exact_mut(n_bytes).nth(idx).unwrap();

        Some(PlacementMut::new(bytes, n_beads))
    }

    pub fn get(&self, pos: UVec3) -> Option<Placement> {
        if pos.x >= self.dimensions.x || pos.y >= self.dimensions.y || pos.z >= self.dimensions.z {
            return None;
        }

        // TODO: Is it wasteful to 'compute' this here? Or does it reduce down well?
        let n_beads = self.solvent.beads.len();
        let n_bytes = n_beads.div_ceil(8);
        let idx = index_3d(pos, self.dimensions);
        // This unwrap should be safe, since we checked at the start of the function.
        let bytes = self.placements.chunks_exact(n_bytes).nth(idx).unwrap();

        Some(Placement::new(bytes, n_beads))
    }
}

pub struct Placement<'ps> {
    bytes: &'ps [u8],
    len: usize,
    idx: usize,
}

impl<'ps> Placement<'ps> {
    fn new(bytes: &'ps [u8], len: usize) -> Self {
        Self { bytes, len, idx: 0 }
    }
}

impl Iterator for Placement<'_> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.len {
            return None;
        }

        let byte = self.bytes[self.idx / 8];
        let bit = (byte >> (self.idx % 8)) & 1;
        self.idx += 1;
        Some(bit > 0)
    }
}

pub struct PlacementMut<'ps> {
    bytes: &'ps mut [u8],
    len: usize,
}

impl<'ps> PlacementMut<'ps> {
    fn new(bytes: &'ps mut [u8], len: usize) -> Self {
        Self { bytes, len }
    }

    /// Set the bit value for a solvent bead to represent that its spot is occupied.
    ///
    /// This communicates that a solvent particle cannot be rendered at this spot.
    pub fn set_occupied(&mut self, idx: usize) {
        if idx >= self.len {
            panic!(
                "index out of range (length was {} but index was {idx})",
                self.len
            )
        }

        let byte = &mut self.bytes[idx / 8];
        *byte |= 1 << (idx % 8);
    }
}

pub struct Cookies {
    dimensions: UVec3,
    cookie_size: Vec3,
    cookies: Box<[Box<[Vec3]>]>,
}

impl Cookies {
    pub fn new(structure: &Structure, cookie_size: Vec3, dimensions: UVec3) -> Self {
        let n_cells = dimensions.x as usize * dimensions.y as usize * dimensions.z as usize;
        let mut cookies = vec![Vec::new(); n_cells];

        // TODO: Great target for parallelization! (If at all necessary?)
        // Place the bead's position in the appropriate cookie.
        // Like raisins or pieces of chocolate. But in our case the cookies aren't round but
        // cuboid. Look, I don't even remember why I picked this name.
        for bead in &structure.beads {
            let pos = bead.pos;
            let cell_pos = (pos / cookie_size).floor().as_uvec3();
            let idx = index_3d(cell_pos, dimensions);
            cookies[idx].push(pos);
        }

        Self {
            dimensions,
            cookie_size,
            cookies: cookies
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
        cell_pos.as_vec3() * self.cookie_size
    }
}
