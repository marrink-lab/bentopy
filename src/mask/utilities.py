from pathlib import Path

import MDAnalysis as mda
import numpy as np


def voxels_to_gro(path: Path, arr, scale: float = 1.0, place_in_center: bool = True):
    scale *= 10.0  # We go from Å to nm scale.
    n_atoms = np.product(arr.shape)

    # Set up our output universe.
    u = mda.Universe.empty(
        n_atoms,
        n_residues=n_atoms,
        atom_resindex=range(n_atoms),
        residue_segindex=[0] * n_atoms,
        trajectory=True,
    )

    # Fill it up with an 3×n_atoms coordinate array.
    indices = np.array(list(np.ndindex(arr.shape)))
    u.atoms.positions = indices * scale

    if place_in_center:
        # Place the atoms in the _middle_ of their voxel, rather than at the voxel origin.
        u.atoms.positions += scale / 2

    # Derive the box size from the array shape.
    box = np.array(arr.shape) * scale
    u.dimensions = np.array([*box, 90, 90, 90], dtype=np.int32)

    # Assign appropriate residue names. We want to distinguish between the
    # different values we find in the array (the different compartments!).
    names = [str(name) for name in arr.flatten()]
    u.add_TopologyAttr("name", names)

    # Write them!
    u.atoms.write(path)
