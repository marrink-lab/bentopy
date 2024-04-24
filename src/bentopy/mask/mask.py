from pathlib import Path

import MDAnalysis as mda
import numpy as np
from mdvcontainment import Containers

# Import .config.setup_parser here so it is accessible through mask.setup_parser.
from .config import setup_parser as setup_parser
from .utilities import voxels_to_gro

log = print


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

    # Write the plot of the containments if desired.
    if args.plot is None and args.interactive:
        log("Do you want to write an interactive plot to inspect the containments?")
        log("Leave empty to skip or provide an output path.")
        while True:
            path = input("(html) -> ").strip()
            if len(path) == 0:
                plot_path = None
                break
            plot_path = Path(path)
            if plot_path.suffix == ".html":
                break
            log(f"Path must have an html extension. Found '{plot_path.suffix}'.")
    else:
        plot_path = args.plot
    if plot_path is not None:
        log(f"Plotting to {plot_path}.")
        containment.plot(name=str(plot_path))

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

    mask_res = args.mask_resolution
    mask_shape = zoomed.shape
    mask_size = tuple(np.array(mask_shape) / mask_res)
    log(f"Size of final voxel mask is {mask_shape} at a {mask_res} nm resolution.")
    log(f"This corresponds to a final mask size of {mask_size} nm.")
    if args.debug_voxels is not None:
        voxels_path = args.debug_voxels
        log(f"Writing mask voxels debug file to {voxels_path}... ", end="")
        voxels_to_gro(voxels_path, zoomed)
        log("done.")
    log(f"Writing the voxel mask to {args.output}... ", end="")
    np.savez(args.output, zoomed)
    log("done.")
