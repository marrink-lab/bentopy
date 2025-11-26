use std::{cell::RefCell, collections::HashMap};

use bentopy::core::config::CompartmentID;

use crate::mask::{Mask, distance_mask_grow};

type DistanceMasksKey = u64;

pub struct Compartment {
    pub id: CompartmentID,
    pub mask: Mask,
    // We make use of some internal mutability, here, through the RefCell. This allows us to lazily
    // set up the distance masks, such that they are created only once we need them.
    pub(crate) distance_masks: RefCell<HashMap<DistanceMasksKey, Mask>>,
}

impl Compartment {
    /// Get or create and return a cloned distance mask for some voxel distance.
    ///
    /// Note that the distance is in terms of voxels, not in nm.
    pub fn get_distance_mask(&self, distance: f32) -> Mask {
        let key = distance as u64;
        {
            // In this scope, we mutably borrow the distance_masks. At the end of the scope, the
            // mutable borrow is dropped, and we are okay to perform a non-mutable borrow after it.
            let mut distance_masks = self.distance_masks.borrow_mut();
            distance_masks
                .entry(key)
                .or_insert_with(|| distance_mask_grow(&self.mask, distance as u64));
        }
        // We can safely get at the `key` because if it was empty, we have inserted it above.
        // We clone here because we can't guarantee that the value is not moved by a resize of
        // the HashMap.
        self.distance_masks.borrow().get(&key).unwrap().clone()
    }
}
