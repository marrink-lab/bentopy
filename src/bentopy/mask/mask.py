from pathlib import Path

import MDAnalysis as mda
import numpy as np
from mdvcontainment import Containers

# Import .config.setup_parser here so it is accessible through mask.setup_parser.
from .config import setup_parser as setup_parser
from .utilities import voxels_to_gro

log = print


def main(args):
    # First, we want to see whether everything is set in case --no-interactive is used.
    if not args.interactive:
        log("Running in non-interactive mode.")
        if args.labels is None:
            log("ERROR: No labels were specified.", end=" ")
            log("In non-interactive mode, at least one label must be provided.")
            exit(1)

    # Read in structures.
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

    # Get the label array that forms the basis of our masking.
    label_array = containment.data["relabeled_combined_label_array"]
    # Write the labels to a gro file if desired.
    if args.inspect_labels_path is None and args.interactive:
        log("Do you want to write a label map as a gro file to view the containments?")
        log("Leave empty to skip or provide an output path.")
        while True:
            path = input("(gro) -> ").strip()
            if len(path) == 0:
                labels_path = None
                break
            labels_path = Path(path)
            if labels_path.suffix == ".gro":
                break
            log(f"Path must have a gro extension. Found '{labels_path.suffix}'.")
    else:
        labels_path = args.inspect_labels_path
    if labels_path is not None:
        log(f"Writing labels voxels debug file to {labels_path}... ", end="")
        voxels_to_gro(labels_path, label_array)
        log("done.")

    # Let's select our labels.
    possible_labels = np.unique(label_array)
    if args.interactive and args.labels is None:
        log("No compartment labels have been selected, yet.", end=" ")
        log("Select one or more to continue.")
        log(f"Options: {possible_labels}.", end=" ")
        log("Provide them as a space-separated list followed by a return.")
        labels = []
        while len(labels) == 0:
            try:
                labels.extend(map(lambda label: int(label), input("-> ").split()))
            except ValueError:
                log("Could not parse the labels. Make sure they are well-formed.")
            for label in labels:
                if label not in possible_labels:
                    log(f"'{label}' is not a valid compartment label.", end=" ")
                    log("Please try again.")
                    labels.clear()
                    break
    else:
        # Making sure we this is the case, even though we do this check at the top as well.
        assert (
            args.labels is not None
        ), "No labels are specified and the mode is non-interactive so we can't ask."
        labels = args.labels
        for label in labels:
            if label not in possible_labels:
                log(f"ERROR: '{label}' is not a valid compartment label.")
                exit(1)
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
