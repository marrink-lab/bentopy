import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from .utilities import voxels_to_gro
from mdvcontainment import Containers

DEFAULT_CONTAINMENT_RESOLUTION = 1.0
DEFAULT_MASK_RESOLUTION = 0.5

log = print


def process_labels(labels: str) -> list[int]:
    return [int(label) for label in labels.split()]


def setup_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="Set up masks based on structures and compartment segmentations.",
            prog="mask",
            epilog='"Relabel, relabel, relabel."',
        )
    parser.add_argument(
        "input",
        type=Path,
        help="Input structure file to subject to segmentation.",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Output path for the resulting voxel mask.",
    )
    parser.add_argument(
        "--no-interactive",
        dest="interactive",
        action="store_false",
        help="""Do not ask for command line input at runtime. May be desirable 
        when all inputs are known and in a scripted context.
        (normally interactive)""",
    )
    parser.add_argument(
        "--containment-resolution",
        default=DEFAULT_CONTAINMENT_RESOLUTION,
        type=float,
        help="Resolution for the compartment finding routine (nm). (default: %(default)s nm)",
    )
    parser.add_argument(
        "--mask-resolution",
        default=DEFAULT_MASK_RESOLUTION,
        type=float,
        help="Voxel size (resolution) for the exported mask (nm). (default: %(default)s nm)",
    )
    parser.add_argument(
        "--selection",
        default="not resname W ION",
        type=str,
        help="MDAnalysis selection string for the atom group over which to perform the segmentation. (default: '%(default)s')",
    )
    parser.add_argument(
        "--labels",
        type=process_labels,
        help="Pre-selected compartment labels, comma-separated.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        nargs="?",
        const="containment.html",
        help="""Optional output path for the visualization page generated by 
        mdvcontainment during segmentation. 
        (when used, default: %(const)s).""",
    )
    parser.add_argument(
        "-b",
        "--inspect-labels-path",
        type=Path,
        nargs="?",
        const="labels.gro",
        help="""Optional output path to a gro file to write labeled voxel positions to.
        This file can be inspected with molecule viewers, which is very helpful 
        in determining which labels match the compartments you want to select.
        (when used, default: %(const)s)""",
    )
    parser.add_argument(
        "--debug-voxels",
        type=Path,
        help="""Write the final voxel mask as a gro file for inspection with molecule viewers. 
        This can be useful when you want to verify the voxel mask that is produced for some selection of labels.""",
    )
    return parser


def main(args):
    structure_path = args.input
    log(f"Reading in structure from {structure_path}... ", end="")
    u = mda.Universe(structure_path)
    log(f"done. (Read {u.atoms.n_atoms} atoms.)")

    selection = u.select_atoms(args.selection)
    log(f"Selected {selection.n_atoms} atoms according to '{args.selection}'.")

    log("Calculating containment...")
    log("\n--- mdvcontainment --- 8< ---")
    containment = Containers(selection, args.containment_resolution)
    log("--- >8 --- mdvcontainment ---\n")

    log("Found the following components:")
    log(f"\troot:\t{containment.get_root_components()}")
    log(f"\tleaf:\t{containment.get_leaf_components()}")
    # TODO: Add an interactive question about whether to plot and if so what path if it is not provided through the args.
    if args.plot is not None:
        log(f"Plotting to {args.plot}.")
        containment.plot(name=str(args.plot))

    label_array = containment.data["relabeled_combined_label_array"]
    if args.inspect_labels_path is not None:
        labels_path = args.inspect_labels_path
        log(f"Writing labels voxels debug file to {labels_path}... ", end="")
        voxels_to_gro(labels_path, label_array)
        log("done.")

    # TODO: Is there a neater way of getting the different labels from the Containment?
    possible_labels = np.unique(label_array)

    # TODO: Let the behavior depend on whether --interactive is set. Now we will always ask for input().
    if args.labels is None:
        log(
            "No compartment labels have been selected, yet. Select one or more to continue."
        )
        log(
            f"Options: {possible_labels}. Provide them as a space-separated list followed by a return."
        )
        labels = []
        while len(labels) == 0:
            labels.extend(map(lambda label: int(label), input("-> ").split()))
            for label in labels:
                if label not in possible_labels:
                    log(
                        f"'{label}' is not a valid compartment label. Please try again."
                    )
                    labels.clear()
                    break
    else:
        labels = args.labels
    log(f"Selected the following labels: {labels}.")

    # TODO: Reconsider if args.output should be mandatory when --interactive is also there.
    if args.output is None:
        log("Please specify an output path or name for the mask (npz).")
        mask_path = ""
        while len(mask_path) == 0:
            mask_path = input("-> ").trim()
    else:
        mask_path = args.output

    compartment = np.isin(label_array, labels)
    # TODO: Display this as % voxels as placeable space.
    (_things, (full, free)) = np.unique(compartment, return_counts=True)
    containment_voxel_volume = args.containment_resolution**3
    log("Selected compartment contains:")
    log(f"    occupied:\t{full} voxels\t({full * containment_voxel_volume:.1f} nm³)")
    log(f"   available:\t{free} voxels\t({free * containment_voxel_volume:.1f} nm³)")

    # HACK: We are just assuming that args.containment_resolution and args.mask_resolution are integer-divisible by each other.
    zoom = int(args.containment_resolution / args.mask_resolution)
    log(f"Output mask resolution is set to {args.mask_resolution} nm.")
    log(f"Zoom factor from containment voxels to mask voxels is {zoom}.")
    zoomed = compartment.repeat(zoom, axis=0).repeat(zoom, axis=1).repeat(zoom, axis=2)

    log(f"Size of final voxel mask is {zoomed.shape}.")
    if args.debug_voxels is not None:
        voxels_path = args.debug_voxels
        log(f"Writing mask voxels debug file to {voxels_path}... ", end="")
        voxels_to_gro(voxels_path, zoomed)
        log("done.")
    log(f"Writing the voxel mask to {mask_path}... ", end="")
    np.savez(mask_path, zoomed)
    log("done.")
