import numpy as np

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
        self.compartments = []
        for c in compartments:
            id = c["id"]
            if "shape" in c:
                compartment = Compartment(id, shape=c["shape"])
            elif "voxels" in c:
                compartment = Compartment(id, voxels=c["voxels"])
            self.compartments.append(compartment)

    def background(self, compartment_ids=[], onto=None):
        # Adjust size and padding for resolution.
        size = (np.array(self.size) / self.resolution).astype(int)
        if onto is None:
            background = np.ones(size, dtype=np.float32)
        else:
            background = onto

        width, height, depth = size
        for compartment in filter(lambda c: c.id in compartment_ids, self.compartments):
            mask = compartment.mask(width, height, depth, padding=0)
            # Apply the mask
            background[mask] = 0
        return background
