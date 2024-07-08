import argparse
from pathlib import Path
from time import time
from sys import stderr

import freud
import MDAnalysis as mda
import numpy as np

def log(*args, **kwargs):
    print(*args, **kwargs, file=stderr)


def check(args):
    # Read in structure.
    structure_path = args.file
    log(f"Reading in structure from {structure_path}... ", end="")
    u = mda.Universe(structure_path)
    log(f"done. (Read {u.atoms.n_atoms} atoms.)")

    # Adjust the box for Freud and convert the cutoff distance.
    box = mda_box_to_freud(u.dimensions)
    cutoff_nm = args.cutoff
    cutoff_A = cutoff_nm * 10.0  # Convert from nm to Å.

    # Set up the locality query.
    positions = u.atoms.positions
    aq = freud.locality.AABBQuery(box, positions)
    neighbors = aq.query(positions, {"r_max": cutoff_A}).toNeighborList()
    log(f"The AABB query was set up with a {cutoff_nm} nm cutoff.")

    # Sift through the neighbors to find collisions.
    collision = False
    min_distance = None
    collisions = [] if args.output_collisions else None
    try:
        for id, hit in neighbors:
            a = u.atoms[id]
            b = u.atoms[hit]
            if a.resid != b.resid:
                if args.ignore_same_resname and a.resname == b.resname:
                    continue
                collision = True
                distance = np.linalg.norm(a.position - b.position) / 10.0  # From Å to nm.
                new_low = min_distance is not None and min_distance > distance
                if new_low or min_distance is None:
                    min_distance = distance
                print(
                    f"Distance between atom {id:>6} and {hit:>6} (residue {a.resname} ({a.resid:>3}) and {b.resname} ({b.resid:>3})) is {distance:<6.3} <= {cutoff_nm} nm.",
                    "(new smallest distance)" if new_low else "",
                )
                if collisions is not None:
                    collisions.append((a.position + b.position) / 2)
                if args.exit_early:
                    break
        log("Done.")
    except KeyboardInterrupt:
        log("Stopping the search.")

    # If desired, write out a structure with beads at the collision sites.
    if args.output_collisions is not None:
        log(f"Writing collision coordinates to {args.output_collisions}... ", end="")
        start = time()
        natoms = len(collisions)
        positions = np.array(collisions)
        assert positions.shape == (natoms, 3)
        # Create a new Universe. We need to set trajectory to True in order to add positions.
        cu = mda.Universe.empty(natoms, trajectory=True)
        cu.atoms.positions = positions
        cu.dimensions = u.dimensions
        cu.atoms.write(args.output_collisions)
        duration = time() - start
        log(f"Done in {duration:.3} s")

    # Report whether we found any collisions through the error number.
    return 1 if collision else 0


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


def main():
    parser = argparse.ArgumentParser(
        description="Check for collisions in rendered structures.",
        prog="bentopy-check",
    )
    parser.add_argument("file", type=Path, help="File to check (a structure file).")
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.4,
        help="Collision distance between beads in nm. (default: %(default)s nm)",
    )
    parser.add_argument(
        "--output-collisions",
        type=Path,
        help="Write out a structure file with dummy beads at the collision sites.",
    )
    parser.add_argument(
        "--ignore-same-resname",
        action="store_true",
        help="Ignore collisions between particles with the same residue name.",
    )
    parser.add_argument(
        "-e",
        "--exit-early",
        action="store_true",
        help="Exit early at the first collision.",
    )

    args = parser.parse_args()
    return check(args)


if __name__ == "__main__":
    main()
