import numpy as np
import MDAnalysis as MDA


class Segment:
    def __init__(
        self,
        name,
        target_number,
        path,
        resolution,
        compartment_ids,
        rotation_axes,
        center,
    ):
        self.name = name
        self.target_number = target_number
        self.path = path
        self.rotation = np.eye(3)  # We start out with the identity matrix.
        self.resolution = resolution
        self.compartment_ids = compartment_ids
        if rotation_axes is None:
            self.rotation_axes = parse_rotation_axes("xyz")
        else:
            self.rotation_axes = parse_rotation_axes(rotation_axes)
        if center is None:
            self.center = parse_center("auto, auto, auto")
        else:
            self.center = parse_center(center)
        self._points = None
        # TODO: Perhaps this ought to be a dictionary, but then we would need a
        # method of normalizing rotation matrices to a lower precision in order
        # to make the bit patterns of ~identical rotations exactly the same.
        # Otherwise there'd be no benefit.
        self.batches = []  # A batch is a set of placements with a particular rotation.

    def points(self):
        """
        Returns the centered point cloud of the atoms in this Segment.

        Positions are given in nanometers.

        If no point cloud is stored internally, yet, it will be loaded from
        this Segment's path.
        """
        if self._points is None:
            points = load_points(self.path)
            center = np.mean(points, axis=0)
            # Apply the configured center correction.
            points -= center + self.center_translation()
            self._points = points

        return self._points

    def voxels(self):
        """
        Returns the voxels for this Segment's point cloud, with this
        Segment's rotation applied before voxelization.
        """
        return voxelize(self.rotation @ self.points().T, self.resolution, tighten=True)

    def center_translation(self, resolution=1.0, dtype=float):
        return np.array(
            [0.0 if c is None else c / resolution for c in self.center], dtype=dtype
        )

    def set_rotation(self, rotation):
        """
        Set the internal rotation by providing a rotation matrix.

        The provided rotation matrix will be processed according to the
        internal rotation axes constraints. For example, if the only rotational
        axis that is set is z, only the z component of the rotation will be
        stored.

        If no rotation axes are set, the internal rotation field will remain
        the same.
        """
        # If all axes are accepted for rotations, we can just use the new rotation directly.
        if all(self.rotation_axes):
            self.rotation = rotation
            return

        # Otherwise, deconstruct the rotation and build it up according to the set axes.
        from scipy.spatial.transform import Rotation as R

        r = R.from_matrix(rotation)
        axis_rotations = r.as_euler("xyz")  # Decompose into axis angles.
        # Finally, we zero out the ommitted rotation axes by multiplying these
        # axis rotations with the rotation axes booleans.
        axis_rotations_constrained = axis_rotations * self.rotation_axes
        self.rotation = R.from_euler("xyz", axis_rotations_constrained).as_matrix()

    def add_rotation(self, placements):
        """
        Add a batch of placements with a particular rotation.

        The rotation is provided as an `[x, y, z]` list where the first value
        represents the rotation in degrees around the _x_-axis, etc.
        """
        self.batches.append((self.rotation.tolist(), placements))


def parse_center(center: str):
    """
    Parse a center configuration string.

    For example:
    - `'auto, auto, auto' -> (None, None, None)`
    - `'auto, auto, 5.1' -> (None, None, 5.1)`
    - `'1.2, 3.4, 5.6' -> (1.2, 3.4, 5.6)`
    """
    axes = [a.strip() for a in center.split(",")]
    return tuple(None if a == "auto" else float(a) for a in axes)


def parse_rotation_axes(axes: str) -> np.ndarray:
    return np.array(("x" in axes, "y" in axes, "z" in axes))


def load_points(path):
    """
    Load atoms from a molecule file at `path` and return them as a point cloud.

    Positions in nanometers.
    """
    # TODO: Revisit the fact that we're creating the universe two times over
    # the runtime (depending on whether the placement list is rendered to a gro
    # file internally).
    # TODO: Consider whether it should be loaded into memory explicitly.
    u = MDA.Universe(path)
    # Convert from Å to nm.
    return u.atoms.positions / 10.0


def voxelize(points, resolution, tighten=False):
    """
    Return the voxels corresponding to a point cloud.

    The number of points must be at least 1. Since there are no sensible return
    values for that case, we throw a ValueError.

    Any voxels returned by this function are quaranteed to have at least one
    filled-in voxel.

    The value of the provided resolution will be the final voxel size.
    """
    _, npoints = points.shape
    if npoints == 0:
        raise ValueError("Cannot meaningfully voxelize a structure without any points.")
    # From here on out, we know that there is at least one point. We will use that knowledge to check our work.

    points = points / resolution
    mins = points.min(axis=1)
    points -= mins[:, None]

    maxs = np.ceil(points.max(axis=1)).astype(int)
    voxels = np.zeros(tuple(max(1, m) for m in maxs))  # FIXME: Better way?
    assert all(
        v >= 1 for v in voxels.shape
    ), f"The voxels array should have a shape of at least (1, 1, 1). It was {voxels.shape}"

    # TODO: This seems silly. Must be a better way.
    for x, y, z in points.T:
        voxels[int(x), int(y), int(z)] = 1

    if tighten:
        # And now, "tighten up those lines!"
        # TODO: We can do this in a more pretty manner, I'm sure.
        mins = np.min(np.where(voxels), axis=1)
        voxels = voxels[mins[0] :, mins[1] :, mins[2] :]
        maxs = np.max(np.where(voxels), axis=1) + 1
        voxels = voxels[: maxs[0], : maxs[1], : maxs[2]]

    assert (
        np.sum(voxels) > 0
    ), f"There is at least one point in this structure, yet not a single voxel was filled in."

    return voxels
