import warnings
import pickle
from time import time
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from mdvcontainment import Containment

from .config import setup_parser
from .utilities import voxels_to_gro

# Let's ignore the wordy warnings we tend to get from MDAnalysis.
warnings.filterwarnings("ignore")

log = print


def mask(args):
    # First, we want to verify our arguments.
    # Down the line, we assume that the containment resolution can be treated
    # as integer multiple of the mask resolution. Let's check that at the top.
    if args.containment_resolution % args.mask_resolution > 0.01:
        log("ERROR: The containment and mask resolutions are not well-formed.")
        log(
            f"The containment resolution ({args.containment_resolution})",
            f"must be a multiple of the mask resolution ({args.mask_resolution}),",
            "such that `containment_resolution % mask resolution == 0`.",
        )
        return 1
    if args.mask_resolution < 0:
        log(f"ERROR: The mask resolution ({args.mask_resolution}) cannot be negative.")
        return 1
    zoom = int(args.containment_resolution / args.mask_resolution)
    # Everything must be set correctly in case --no-interactive is used.
    if not args.interactive:
        log("Running in non-interactive mode.")
        if args.output is None:
            log("WARNING: No mask output path was specified.")
            log("The computed mask will not be written to disk.")
            log("Like tears in rain.")
        if args.labels is None and not args.autofill:
            log("ERROR: No labels were specified.")
            log("In non-interactive mode, at least one label must be provided manually (`--label`) or automatically (`--autofill`).")
            return 1

    # Read in structures from a structure file or from a cache..
    u = None
    caching = not args.no_cache  # Whether caching is enabled.
    structure_path = args.input
    if caching:
        cached_path = resolve_cache_path(structure_path)

        # Try to find a cached equivalent to the structure.
        cache_exists = cached_path.is_file()
        if cache_exists:
            log(f"Reading in cached structure from {cached_path}...")
            log("You can reload a changed structure by removing the cached file.")
            log("Caching can be disabled with --no-cache.")
            with open(cached_path, "rb") as cache_file:
                start = time()
                u = pickle.load(cache_file)
                dur = time() - start
                log(
                    f"Done reading the cached file. (Read {u.atoms.n_atoms} atoms in {dur:.1f} s.)"
                )

    # If the file was not cached, we read it from the structure file.
    if u is None:
        # Read the structure file.
        log(f"Reading in structure from {structure_path}... ", end="")
        start = time()
        u = mda.Universe(structure_path)
        dur = time() - start
        log(f"done. (Read {u.atoms.n_atoms} atoms in {dur:.1f} s.)")

        # Cache the universe if desired and necessary.
        if caching and not cache_exists:
            with open(cached_path, "wb") as cache_file:
                pickle.dump(u, cache_file)

    # Apply the selection before we hand it to mdvcontainment.
    selection = u.select_atoms(args.selection)
    log(f"Selected {selection.n_atoms} atoms according to '{args.selection}'.")

    # Calculate the containments. This is all in the hands of mdvcontainment.
    log("Calculating containment... ", end="")
    start = time()
    containment = Containment(
        selection,
        resolution=args.containment_resolution,
        closure=args.closing,
        slab=args.slab,
        verbose=args.verbose,
        no_mapping=True, # Mapping takes some time and is not used at all in this context.
        betafactors=False, # Betafactors are not used in this context and are pretty slow to instantiate.
        
    )
    duration = time() - start
    log(f"Done in {duration:.3} s.")

    # Show what we found.
    log("Found the following node groups:")
    log(f"        root:\t{npc(containment.voxel_containment.root_nodes)}")
    log(f"        leaf:\t{npc(containment.voxel_containment.leaf_nodes)}")
    # And show it as a tree of containments.
    # TODO: This has got to change on Bart's end. I just want a string that I can print at my own leisure. @Bart
    containment.voxel_containment.print_containment()

    # Write the plot of the containments if desired.
    show_plot = args.plot
    if show_plot is None and args.interactive:
        log("Do you want to inspect a plot of the containments?")
        while True:
            answer = input("[y/N] -> ").strip()
            if len(answer) == 0 or answer.lower() == "n":
                break
            if answer == "y":
                show_plot = True
                break
    elif isinstance(args.plot, Path):
        # TODO: Write the plot to a file. @Bart
        log("WARNING: This is not implemented, yet.")
        log(f"Writing the containment plot to {args.plot}.")
        containment.voxel_containment.draw(file=args.plot)
    if show_plot:
        log("Showing the containment plot.")
        containment.voxel_containment.draw()

    # Get the label array.
    label_array = containment.voxel_containment.components_grid
    # Write the labels to a gro file if desired.
    if args.inspect_labels_path is None and args.interactive:
        log("Do you want to write a label map as a gro file to view the containments?")
        log("Provide an output path. To skip this step, leave this field empty.")
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
    if args.interactive and args.labels is None and args.autofill is False:
        log(
            "No compartment labels have been selected, yet.",
            "Select one or more to continue.\n",
            f"Options: {possible_labels}.",
            "Provide them as a space-separated list followed by a return.",
        )
        labels = []
        while len(labels) == 0:
            try:
                labels.extend(map(lambda label: int(label), input("-> ").split()))
            except ValueError:
                log("Could not parse the labels. Make sure they are well-formed.")
            for label in labels:
                if label not in possible_labels:
                    log(f"'{label}' is not a valid compartment label.")
                    log("Please try again.")
                    labels.clear()
                    break
    elif args.autofill is True:
        labels = containment.voxel_containment.leaf_nodes
        log(f"Automatically choosing leaf components: {npc(labels)}.")
    else:
        # Making sure we this is the case, even though we do this check at the top as well.
        assert (
            args.labels is not None
        ), "No labels are specified and the mode is non-interactive so we can't ask."
        labels = args.labels
        for label in labels:
            if label not in possible_labels:
                log(f"ERROR: '{label}' is not a valid compartment label.")
                return 1
    log(f"Selected the following labels: {npc(labels)}.")

    # Get our compartment by masking out all voxels that have our selected labels.
    compartment = containment.voxel_containment.get_voxel_mask(labels)
    full = np.count_nonzero(compartment == True)
    free = np.count_nonzero(compartment == False)
    containment_voxel_volume = args.containment_resolution**3
    full_volume = full * containment_voxel_volume
    free_volume = free * containment_voxel_volume
    total_volume = full_volume + free_volume
    full_frac = full_volume / total_volume
    free_frac = free_volume / total_volume
    log("Selected compartment contains:")
    log(f"    occupied:\t{full} voxels\t({full_frac:.1%}, {full_volume:.1f} nm³)")
    log(f"   available:\t{free} voxels\t({free_frac:.1%}, {free_volume:.1f} nm³)")

    # Produce our final output mask according to the specified output mask resolution.
    log(f"Output mask resolution is set to {args.mask_resolution} nm.")
    log(f"Zoom factor from containment voxels to mask voxels is {zoom}.")
    zoomed = compartment.repeat(zoom, axis=0).repeat(zoom, axis=1).repeat(zoom, axis=2)

    # Report a summary of the final mask's dimensions.
    mask_res = args.mask_resolution
    mask_shape = zoomed.shape
    mask_size = tuple(npc(np.array(mask_shape) * mask_res))
    log(f"Size of final voxel mask is {mask_shape} at a {mask_res} nm resolution.")
    log(f"This corresponds to a final mask size of {mask_size} nm.")

    # Write out a debug voxels gro of the mask if desired.
    if args.debug_voxels is None and args.interactive:
        log("Do you want to write the voxel mask as a gro file to inspect it?")
        log("Warning: This file may be quite large, depending on the mask resolution.")
        log("Provide an output path. To skip this step, leave this field empty.")
        while True:
            path = input("(gro) -> ").strip()
            if len(path) == 0:
                voxels_path = None
                break
            voxels_path = Path(path)
            if voxels_path.suffix == ".gro":
                break
            log(f"Path must have a gro extension. Found '{voxels_path.suffix}'.")
    else:
        voxels_path = args.debug_voxels
    if voxels_path is not None:
        log(f"Writing mask voxels debug file to {voxels_path}... ", end="")
        voxels_to_gro(voxels_path, zoomed)
        log("done.")

    # Determine the voxel mask output path.
    output_path = args.output
    if args.output is None and args.interactive:
        log("Please provide an output path for the final voxel mask.")
        while True:
            path = input("(npz) -> ").strip()
            if len(path) == 0:
                output_path = None
                break
            output_path = Path(path)
            if output_path.suffix == ".npz":
                break
            log(f"Path must have an npz extension. Found '{output_path.suffix}'.")

    # Finally, write out the voxel mask when desired.
    if output_path is not None:
        log(f"Writing the voxel mask to {output_path}... ", end="")
        np.savez(output_path, zoomed)
        log("done.")


def resolve_cache_path(structure_path):
    """
    Determine what the cache path for the provided structure path is.
    """
    structure_name = structure_path.name
    cached_name = f".cached_{structure_name}.pickle"
    return structure_path.with_name(cached_name)


def npc_single(n):
    match type(n):
        case np.float32 | np.float64:
            return float(n)
        case np.int32 | np.int64:
            return int(n)


def npc(ns):
    "Canonicalize numpy number types to print them properly."
    return [npc_single(n) for n in ns]


def main():
    parser = setup_parser()
    args = parser.parse_args()
    return mask(args)


if __name__ == "__main__":
    main()
