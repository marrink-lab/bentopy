use glam::U64Vec3;
use rand::seq::SliceRandom;

use crate::mask::{Dimensions, Position};
use crate::state::{Compartment, Rng, Space, Voxels};
use crate::Location;

pub struct Session<'s> {
    inner: &'s mut Space,
    locations: &'s mut Locations,
    target: usize,
}

impl<'s> Session<'s> {
    pub fn new(inner: &'s mut Space, locations: &'s mut Locations, target: usize) -> Self {
        Self {
            inner,
            locations,
            target,
        }
    }
}

impl Session<'_> {
    pub fn resolution(&self) -> f32 {
        self.inner.resolution
    }

    pub fn dimensions(&self) -> Dimensions {
        self.inner.dimensions
    }

    pub fn compartments(&self) -> &[Compartment] {
        &self.inner.compartments
    }

    pub fn periodic(&self) -> bool {
        self.inner.periodic
    }

    /// Returns a location if available.
    ///
    /// This function also takes a reference to a random number generator, since the internal
    /// [`Locations`] type may shuffle its indices on demand.
    ///
    /// If the background for this [`Session`] does not have any free locations left, the session
    /// must be ended. This case is indicated with a `None` return value.
    pub fn pop_location(&mut self, rng: &mut Rng) -> Option<Location> {
        // TODO: Consider having the shuffle guess that is passed to pop be dynamic with the number
        // of placements that have already occurred. We could update that number as we go.
        if let Some(location) = self.locations.pop(rng, self.target) {
            return Some(location);
        }

        // We're out of locations. We'll try and renew once.
        self.locations.renew(
            self.inner
                .session_background
                .linear_indices_where::<false>(),
        );

        // If the result is still `None`, there is nothing we can do anymore.
        self.locations.pop(rng, self.target)
    }

    /// Returns `true` if no collisions are encountered between the [`Space`] and the provided
    /// [`Voxels`] at some `position`.
    ///
    /// If a collision is found, the function exits early with a `false` value.
    pub fn check_collisions(&self, voxels: &Voxels, position: Position) -> bool {
        let occupied_indices = voxels.indices_where::<true>().map(U64Vec3::from_array);
        let position = U64Vec3::from_array(position);

        occupied_indices
            .map(|p| (p + position).to_array())
            .all(|idx| !self.inner.session_background.get_periodic(idx))
    }

    pub fn stamp(&mut self, voxels: &Voxels, location: Position) {
        let location = U64Vec3::from_array(location);
        for idx in voxels
            .indices_where::<true>()
            .map(|idx| U64Vec3::from_array(idx) + location)
        {
            self.inner
                .global_background
                .set_periodic(idx.to_array(), true);
            self.inner
                .session_background
                .set_periodic(idx.to_array(), true);
        }
    }

    pub fn exit_session(&mut self) {
        // Not really doing anything here anymore.
    }

    pub const fn position(&self, location: Location) -> Option<Position> {
        self.inner.session_background.spatial_idx(location)
    }
}

impl<'s> Drop for Session<'s> {
    fn drop(&mut self) {
        self.exit_session();
    }
}

pub struct Locations {
    /// Linear indices into the background map.
    indices: Vec<Location>,
    used: usize,
    threshold: usize,

    /// The number of shuffled elements that may be popped from the end of `indices`.
    cursor: usize,
}

impl Locations {
    /// Threshold of spare capacity at which the `locations` [`Vec`] ought to be shrunk.
    const SHRINK_THRESHOLD: usize = (100 * 1024_usize.pow(2)) / std::mem::size_of::<Location>();
    /// Value that determines the threshold at which the locations will be renewed.
    const RENEW_THRESHOLD: f64 = 0.1;

    pub fn new() -> Self {
        Self {
            indices: Vec::new(),
            used: 0,
            threshold: 0,
            cursor: 0,
        }
    }

    /// Renew the locations from an iterator of locations.
    pub fn renew(&mut self, locations: impl Iterator<Item = Location>) {
        let indices = &mut self.indices;
        indices.clear();
        indices.extend(locations);
        // Shrink the `locations` if the spare capacity is getting out of hand.
        if indices.spare_capacity_mut().len() > Self::SHRINK_THRESHOLD {
            indices.shrink_to_fit();
        }
        self.used = 0;
        self.cursor = 0;
        // We want to make sure that the threshold is rounded down, such that it eventually could
        // reach 0, indicating the enclosing Session is exhausted. See `pop`.
        self.threshold = (indices.len() as f64 * Self::RENEW_THRESHOLD).floor() as usize
    }

    fn shuffle(&mut self, rng: &mut Rng, shuffle_guess: usize) {
        let (shuffled, _unshuffled) = self.indices.partial_shuffle(rng, shuffle_guess);
        self.cursor = shuffled.len();
    }

    /// Pop a location.
    ///
    /// This function also takes a reference to a random number generator, since this type may
    /// shuffle its indices on demand. The `shuffle_guess` suggests the number of items that ought
    /// to be shuffled in that scenario.
    ///
    /// A `None` value indicates that this [`Locations`] object must be [renewed](`Self::renew`),
    /// not necessarily that all locations have been popped.
    pub fn pop(&mut self, rng: &mut Rng, shuffle_guess: usize) -> Option<Location> {
        if self.used >= self.threshold {
            return None; // Indicate a renewal is needed.
        }
        if self.cursor == 0 {
            self.shuffle(rng, shuffle_guess);
        }
        let location = self.indices.pop()?;
        self.used += 1;
        assert!(self.cursor > 0);
        self.cursor -= 1;
        Some(location)
    }
}
