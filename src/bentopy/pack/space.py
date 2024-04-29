import numpy as np
from scipy.signal import fftconvolve

from .mask import cuboid_mask, sphere_mask


class Voxels:
    def __init__(self, path):
        self.path = path
        self.voxels_cache = None
        self.mask_cache = None

    def voxels(self):
        if self.voxels_cache is None:
            if self.path.endswith(".npz"):
                with np.load(self.path) as npz:
                    self.voxels_cache = npz[npz.files[0]]  # FIXME: Seems silly.
            else:
                self.voxels_cache = np.load(self.path)
        return self.voxels_cache

    def mask(self, width, height, depth, _padding):
        if self.mask_cache is None:
            voxels = self.voxels()
            space_size = (width, height, depth)
            if voxels.shape != space_size:
                raise ValueError(
                    f"The provided voxel mask has the size {voxels.shape}, but must have the same size as its encapsulating Space {space_size}."
                )
            # TODO: Apply the width, height, depth, padding parameters?
            self.mask_cache = np.where(voxels)
        return self.mask_cache


class Shape:
    kind_spherical = "spherical"
    kind_cuboid = "cuboid"
    kind_none = "none"
    available = [kind_spherical, kind_cuboid, kind_none]

    def __init__(self, kind):
        if kind not in self.available:
            raise ValueError(f"Unknown Shape kind: '{kind}' not in {self.available}")
        self.kind = kind
        self.mask_cache = None

    def mask(self, width, height, depth, padding):
        if self.mask_cache is None:
            if self.kind == self.kind_spherical:
                size = min(width, height, depth)
                # TODO: Make this elliptical or something idk.
                self.mask_cache = np.where(sphere_mask(size, size / 2 - padding))
            elif self.kind == self.kind_cuboid:
                self.mask_cache = cuboid_mask(width, height, depth, padding)
            elif self.kind == self.kind_none:
                self.mask_cache = tuple()
        return self.mask_cache


class Compartment:
    def __init__(self, id, shape=None, voxels=None):
        # Take care of the two incorrect input scenarios.
        both_none = shape is None and voxels is None
        both_some = shape is not None and voxels is not None
        if both_none or both_some:
            raise ValueError("Provide either a `shape` or `voxels` definition.")
        self.id = id
        if shape is not None:
            self.definition = Shape(shape)
        if voxels is not None:
            self.definition = Voxels(voxels["path"])

    def mask(self, width, height, depth, padding):
        return self.definition.mask(width, height, depth, padding)


class Space:
    def __init__(self, size, resolution, compartments):
        assert len(size) == 3, "The size of a space must be 3-dimensional"
        # Size in nm (unaffected by whatever value resolution has).
        self.size = size
        self.resolution = resolution
        # The size, adjusted for the resolution.
        self.effective_size = (np.array(self.size) / self.resolution).astype(int)
        self.compartments = []
        for c in compartments:
            id = c["id"]
            if "shape" in c:
                compartment = Compartment(id, shape=c["shape"])
            elif "voxels" in c:
                compartment = Compartment(id, voxels=c["voxels"])
            self.compartments.append(compartment)

        # The global background remains untouched by any masks and tracks all placements.
        self.global_background = np.zeros(self.effective_size, dtype=np.float32)
        # The session background adopts the relevant masks as well as the placements for some segment.
        self.session_background = None

    def enter_session(self, compartment_ids=[]):
        if self.session_background is not None:
            raise ValueError(
                "Big whoop. Entering a session while the session background was already set."
            )

        self.session_background = np.ones(self.effective_size, dtype=np.float32)
        width, height, depth = self.effective_size
        for compartment in filter(lambda c: c.id in compartment_ids, self.compartments):
            mask = compartment.mask(width, height, depth, padding=0)
            self.session_background[mask] = 0  # Apply the mask.
        self.session_background[self.global_background == 1.0] = 1.0

        indices = np.where(self.session_background == 0.0)
        mins = np.min(indices, axis=1)
        maxs = np.max(indices, axis=1)
        self.squeezed_global_background = self.global_background[
            mins[0] : maxs[0], mins[1] : maxs[1], mins[2] : maxs[2]
        ]
        self.squeezed_session_background = self.session_background[
            mins[0] : maxs[0], mins[1] : maxs[1], mins[2] : maxs[2]
        ]
        self.squeezed_location_offset = mins

    # TODO: Make something with __entry__ here, so it can be used in a with block?
    def exit_session(self):
        if self.session_background is None:
            raise ValueError(
                "Big whoop. Exited the session while the session background was already unset."
            )
        self.session_background = None

    def collisions(self, segment_voxels):
        """
        Identify collisions between the inner session background and the segment voxels using an fft convolution.

        Returns an array of the valid placement points.
        """
        query = np.flip(segment_voxels)
        collisions = fftconvolve(self.squeezed_session_background, query, mode="valid")

        valid_offset = np.array(query.shape) // 2
        # The valid placement points will have a value of 0. Since the floating
        # point operations leave some small errors laying around, we use a quite
        # generous cutoff.
        valid = np.array(np.where(collisions < 1e-4)) + valid_offset[:, None]
        return valid

    def stamp(self, voxels_indices):
        """
        Stamp a set of voxels onto the session and global backgrounds.

        Stomp on its stamps, if you will.
        """

        self.squeezed_global_background[
            voxels_indices[0, :],
            voxels_indices[1, :],
            voxels_indices[2, :],
        ] = 1.0
        self.squeezed_session_background[
            voxels_indices[0, :],
            voxels_indices[1, :],
            voxels_indices[2, :],
        ] = 1.0
