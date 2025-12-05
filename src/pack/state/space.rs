use std::collections::HashSet;

use bentopy::core::config::{CompartmentID, Quantity};

use crate::mask::{Dimensions, Mask};
use crate::session::{Locations, Session};
use crate::state::{Size, compartment::Compartment};

pub type Compartments = Vec<Compartment>;

pub struct Space {
    pub size: Size,
    pub dimensions: Dimensions,
    pub resolution: f32,
    pub compartments: Compartments,
    pub periodic: bool,

    pub(crate) global_background: Mask,
    pub(crate) session_background: Mask,
    /// The previous session's compartment IDs are used to see if a renewal of locations is
    /// necessary between subsequent segment placements.
    ///
    /// When set to `None`, a renewal of locations is due for the next session, regardless of the
    /// previous session's compartment IDs.
    pub(crate) previous_compartments: Option<HashSet<CompartmentID>>,
}

impl Space {
    pub fn enter_session<'s>(
        &'s mut self,
        compartment_ids: impl IntoIterator<Item = CompartmentID>,
        locations: &'s mut Locations,
        quantity: Quantity,
    ) -> Session<'s> {
        let compartment_ids = HashSet::from_iter(compartment_ids);

        // TODO: Consider caching this volume like we do for the same compartments below.
        //       The volume can just be associated with a set of previous compartments.
        let target = quantity.bake(|| self.volume(&compartment_ids));

        // Set up a new session background if necessary.
        // Otherwise, leave the session background and locations alone. The session background
        // will stay exactly the same, since it was already set up for this set of
        // compartments. The locations are likely still valid.
        let same_previous_compartments = self
            .previous_compartments
            .as_ref()
            .is_some_and(|prev| prev == &compartment_ids);
        if !same_previous_compartments {
            // Clone the global background, which has all structures stamped onto it.
            self.session_background = self.global_background.clone();

            // Apply the compartments to the background.
            if let Some(merge) = self
                .compartments
                .iter()
                .filter(|comp| compartment_ids.contains(&comp.id))
                .map(|comp| comp.mask.clone())
                .reduce(|mut acc, mask| {
                    acc.merge_mask(&mask);
                    acc
                })
            {
                self.session_background.apply_mask(&merge);
            }
            self.previous_compartments = Some(compartment_ids);

            // We must renew the locations as well, based on the newly masked session background.
            locations.renew(self.session_background.linear_indices_where::<false>());
        }

        Session::new(
            self,
            locations,
            target
                .try_into()
                .expect("target cannot not exceed system word size"),
        )
    }

    /// Determine the free voxel volume for the specified compartments.
    fn volume(&self, compartment_ids: &HashSet<CompartmentID>) -> f64 {
        let free_voxels = self
            .compartments
            .iter()
            .filter(|comp| compartment_ids.contains(&comp.id))
            .map(|comp| comp.mask.clone())
            .reduce(|mut acc, mask| {
                acc.merge_mask(&mask);
                acc
            })
            .map(|mask| mask.count::<false>())
            .unwrap_or(0);

        let voxel_volume = (self.resolution as f64).powi(3);
        free_voxels as f64 * voxel_volume
    }
}
