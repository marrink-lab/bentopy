import pickle
import warnings
from pathlib import Path
from time import time

import MDAnalysis as mda
import numpy as np
from mdvcontainment import Containment
from mdvcontainment.voxel_logic import voxels_to_universe

from .config import parse_args
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

    # Read in structures from a structure file or from a cache..
    u = None
    caching = not args.no_cache  # Whether caching is enabled.
    structure_path = args.input
    if caching:
        cached_path = resolve_cache_path(structure_path)

        # We need to check if the requested structure file actually exists.
        # We don't want to end up using a cached file for a structure file that
        # does not exist anymore.
        structure_exists = structure_path.is_file()
        # Try to find a cached equivalent to the structure.
        cache_exists = cached_path.is_file()
        if structure_exists and cache_exists:
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

    # MDVContainment only accepts str, not str | None as of 2.0.0a2.
    morph = args.morph
    if morph == None:
        morph = ""

    # Calculate the containments. This is all in the hands of mdvcontainment.
    log("Calculating containment... ", end="")
    start = time()
    containment = Containment(
        selection,
        resolution=args.containment_resolution,
        morph=morph,
        max_offset=0,  # We accept any result of voxelization.
        verbose=args.verbose,
        no_mapping=True,  # Mapping takes some time and is not used at all in this context.
    )
    duration = time() - start
    log(f"Done in {duration:.3} s.")

    if args.min_size is not None:
        log(f"Applying a {args.min_size} nm³ minimum size view.")
        containment = containment.node_view(min_size=args.min_size)

    # Show what we found.
    log("Found the following node groups:")
    log(f"        root:\t{npc(containment.voxel_containment.root_nodes)}")
    log(f"        leaf:\t{npc(containment.voxel_containment.leaf_nodes)}")
    # And show it as a tree of containments.
    print(containment.voxel_containment)

    # Get the label array.
    label_array = containment.voxel_containment.components_grid

    # Write the labels to a gro file if desired.
    labels_path = args.visualize_labels
    if labels_path is not None:
        if args.exclude_outside:
            negative_root_nodes = set(
                rn for rn in containment.voxel_containment.root_nodes if rn < 0
            )
            nodes = set(containment.voxel_containment.nodes) - negative_root_nodes
        else:
            nodes = containment.voxel_containment.nodes
        labels = containment.voxel_containment.components_grid
        labels_u = voxels_to_universe(labels, nodes=list(nodes), universe=u)
        assert labels_u is not None
        assert labels_u.atoms is not None

        if len(labels_u.atoms) > 0:
            log(f"Writing labels voxels debug file to {labels_path}... ", end="")
            labels_u.atoms.write(labels_path)
            log("done.")
        else:
            log(
                f"Could not write labels to {labels_path}, because there are no voxels to write!"
            )

    # Let's select our labels.
    possible_labels = npc(containment.voxel_containment.nodes)
    if args.labels is not None:
        labels_to_masks = args.labels
        # Substitute None (autofill directives), such that:
        #   [([int | None], str)] -> [([int], str)]
        for labels, mask_path in labels_to_masks:
            final_labels = set()
            for label in labels:
                if label is None:
                    # Here we find the labels for 'autofill'.
                    leaf_nodes = containment.voxel_containment.leaf_nodes
                    final_labels.update(leaf_nodes)
                else:
                    final_labels.add(label)
            labels[:] = list(final_labels)
            for label in labels:
                if label not in possible_labels:
                    log(
                        f"ERROR: {label} for '{mask_path}' is not a "
                        "valid compartment label."
                    )
                    return 1
    else:
        log("No label selection was provided.")
        labels = npc(containment.voxel_containment.leaf_nodes)
        # Make sure that we don't spam the entire screen if there are a lot of
        # leaf node labels.
        nlabels = min(len(labels), 5)
        example_labels = ",".join(str(l) for l in labels[::-1][:nlabels])
        log(f"To write a mask for, e.g., the compartments {example_labels}, use")
        log(f"        -l {example_labels}:mask.npz")
        return 0

    reported = False
    for labels, mask_path in labels_to_masks:
        # Get our compartment by masking out all voxels that have our selected labels.
        compartment = containment.voxel_containment.get_voxel_mask(labels)
        # Produce our final output mask according to the specified output mask resolution.
        zoomed = (
            compartment.repeat(zoom, axis=0).repeat(zoom, axis=1).repeat(zoom, axis=2)
        )

        if not reported:
            # Report a summary of the final masks's dimensions.
            log(f"Output mask resolution is set to {args.mask_resolution} nm.")
            log(f"Zoom factor from containment voxels to mask voxels is {zoom}.")
            mask_res = args.mask_resolution
            mask_shape = zoomed.shape
            mask_size = tuple(npc(np.array(mask_shape) * mask_res))
            log(
                f"Size of final voxel masks is {mask_shape} at a {mask_res} nm resolution."
            )
            log(f"This corresponds to a final mask size of {mask_size} nm.")
            reported = True

        log(f"Creating '{mask_path}' with the following labels: {npc(labels)}.")
        full = np.count_nonzero(compartment == True)
        free = np.count_nonzero(compartment == False)
        full_volume = full * containment.voxel_volume
        free_volume = free * containment.voxel_volume
        total_volume = full_volume + free_volume
        full_frac = full_volume / total_volume
        free_frac = free_volume / total_volume
        log("\tSelected compartment contains:")
        log(f"\t    occupied:\t{full} voxels\t({full_frac:.1%}, {full_volume:.1f} nm³)")
        log(f"\t   available:\t{free} voxels\t({free_frac:.1%}, {free_volume:.1f} nm³)")

        # Write out a debug voxels gro of the mask if desired.
        # TODO: Add ability to give custom path here.
        if args.visualize_masks:
            voxels_path = Path(f"{mask_path}.gro")
            log(f"\tWriting mask voxels debug file to {voxels_path}... ", end="")
            voxels_to_universe(zoomed, universe=u, nodes=[True])
            voxels_to_gro(voxels_path, zoomed, scale=args.mask_resolution)
            log("done.")

        # Finally, write out the voxel mask.
        log(f"\tWriting the voxel mask to {mask_path}... ", end="")
        np.savez(mask_path, zoomed)
        log("done.")


def resolve_cache_path(structure_path):
    """
    Determine what the cache path for the provided structure path is.
    """
    structure_name = structure_path.name
    cached_name = f"#cached_{structure_name}.pickle"
    return structure_path.with_name(cached_name)


def npc_single(n):
    match n:
        case np.float32() | np.float64() | float():
            return float(n)
        case np.int32() | np.int64() | int():
            return int(n)
        case _:
            print(f"WARNING: Unknown type for {n}: {type(n)}")


def npc(ns):
    "Canonicalize numpy number types to print them properly."
    return [npc_single(n) for n in ns]


def main():
    args = parse_args()
    return mask(args)


if __name__ == "__main__":
    main()
