import numpy as np
import MDAnalysis as MDA


class Segment:
    def __init__(self, name, target_number, path, resolution, compartment_ids):
        self.name = name
        self.target_number = target_number
        self.path = path
        self.rotation = np.eye(3)  # We start out with the identity matrix.
        self.resolution = resolution
        self.compartment_ids = compartment_ids
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
            self._points = load_points(self.path)

        return self._points

    def voxels(self):
        """
        Returns the voxels for this Segment's point cloud, with this
        Segment's rotation applied before voxelization.
        """
        return voxelize(self.rotation @ self.points().T, self.resolution, tighten=True)

    def add_rotation(self, placements):
        """
        Add a batch of placements with a particular rotation.

        The rotation is provided as an `[x, y, z]` list where the first value
        represents the rotation in degrees around the _x_-axis, etc.
        """
        self.batches.append((self.rotation.tolist(), placements))


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

    The value of the provided resolution will be the final voxel size.
    """
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

    return voxels
