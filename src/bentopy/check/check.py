from sys import exit

import freud
import MDAnalysis as mda
import numpy as np


log = print


def main(args):
    # Read in structure.
    structure_path = args.file
    log(f"Reading in structure from {structure_path}... ", end="")
    u = mda.Universe(structure_path)
    log(f"done. (Read {u.atoms.n_atoms} atoms.)")

    box = mda_box_to_freud(u.dimensions)
    cutoff_nm = args.cutoff
    cutoff_A = cutoff_nm * 10.0  # Convert from nm to Å.

    positions = u.atoms.positions
    aq = freud.locality.AABBQuery(box, positions)
    neighbors = aq.query(positions, {"r_max": cutoff_A}).toNeighborList()
    log(f"The AABB query was set up with a {cutoff_nm} nm cutoff.")

    collision = False
    min_distance = None
    for id, hit in neighbors:
        a = u.atoms[id]
        b = u.atoms[hit]
        if a.resid != b.resid:
            collision = True
            distance = np.linalg.norm(a.position - b.position) / 10.0 # From Å to nm.
            new_low = min_distance is not None and min_distance > distance
            if new_low or min_distance is None:
                min_distance = distance
            print(
                f"Distance between atom {id:>6} and {hit:>6} (residue {a.resname} ({a.resid:>3}) and {b.resname} ({b.resid:>3})) is {distance:<6.3} <= {cutoff_nm} nm.",
                "(new smallest distance)" if new_low else ""
            )
            if args.exit_early:
                break
    exit(1 if collision else 0)


def mda_box_to_freud(mda_box):
    x, y, z, alpha, beta, gamma = mda_box.astype(np.float64)
    cosa = np.cos(np.pi * alpha / 180)
    cosb = np.cos(np.pi * beta / 180)
    cosg = np.cos(np.pi * gamma / 180)
    sing = np.sin(np.pi * gamma / 180)
    zx = z * cosb
    zy = z * (cosa - cosb * cosg) / sing
    zz = np.sqrt(z**2 - zx**2 - zy**2)
    matrix = np.array([[x, 0, 0], [y * cosg, y * sing, 0], [zx, zy, zz]])
    return matrix.astype(np.float32).T
