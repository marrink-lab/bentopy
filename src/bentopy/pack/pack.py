import argparse
import json
import math
import time
import pathlib

import numpy as np
from scipy.spatial.transform import Rotation as R

from .config import Configuration
from extensions._extensions import py_place as place

VERBOSE = False
RNG = np.random.default_rng()


def log(*args, **kwargs):
    """
    Only print to console if verbose output is desired.
    """
    if VERBOSE:
        print(*args, **kwargs)


def placement_location(valid, selection, segment):
    """
    Returns the valid voxel placement position for the segment.
    """
    segment_center = np.array(segment.shape) // 2
    valid_pos = valid[:, selection]
    return valid_pos - segment_center


# def place_naive(
#     space,
#     segment,
#     max_at_once,
#     max_tries,
# ) -> list | None:
#     """
#     Place a number of segments into the background.
#     """
#
#     # TODO: ROTATE IN THIS FUNCTION.
#
#     t0 = time.time()
#     # Get indices of the open voxels.
#     lx, ly, lz = [min(q, b) for q, b in zip(segment.voxels().shape, space.squeezed_session_background.shape)]
#     valid = np.array(np.where(space.squeezed_session_background[:-lx, :-ly, :-lz] == 0))
#     if valid.size == 0:
#         return None
#     t1 = time.time(); print(f"||||||||| t0->t1  {(t1 - t0)*1000:.3} ms")
#
#     placements = []
#     hits = 0
#     tries = 0  # Number of times segment placement was unsuccesful.
#     previously_selected = set()
#     while hits < max_at_once:
#         t2 = time.time()
#         if tries >= max_tries:
#             # Give up.
#             log("% ! tries is exceeding max_tries")
#             break
#         if len(previously_selected) == len(valid[0]):
#             log("% ! all open positions have been tried")
#             break
#         while True:
#             selection = RNG.integers(0, len(valid[0]))
#             if selection not in previously_selected:
#                 previously_selected.add(selection)
#                 break
#         t3p = time.time(); print(f"||||||||| t2->t3p  {(t3p - t2)*1000:.3} ms")
#
#         # Make sure that this placement does not overlap with another
#         # previously selected location.
#         location = valid[:, selection]
#         prospect = np.where(segment.voxels()) + location[:, None]
#         # Check for collisions at the prospective site.
#         free = not np.any(
#             space.squeezed_session_background[
#                 prospect[0, :], prospect[1, :], prospect[2, :]
#             ]
#         )
#         t3 = time.time(); print(f"||||||||| t3p->t3  {(t3 - t3p)*1000:.3} ms")
#
#         if free:
#             space.stamp(prospect)
#
#             # If we are looking at a squeezed background, we need to correct the offset.
#             # We also apply the segment's center translation.
#             corrected_location = (
#                 location + space.squeezed_location_offset
#             ) * space.resolution
#             placements.append(corrected_location.tolist())
#             hits += 1
#         else:
#             tries += 1
#             log(
#                 f"%       {tries = }/{max_tries},\thits = {hits}/{max_at_once}",
#                 end="\r",
#             )
#         t3a = time.time(); print(f"||||||||| t3->t3a  {(t3a - t3)*1000:.3} ms")
#     t4 = time.time(); print(f"||||||||| t1->t4  {(t4 - t1)*1000:.3} ms")
#     log("\x1b[K", end="")  # Clear the line since we used \r before.
#     log(f"% . placed {hits} segments with {tries}/{max_tries} misses")
#
#     return placements


def place_batched(
    space,
    segment,
    max_at_once,
    max_tries,
    threshold_coefficient,
) -> list | None:
    """
    Place a number of segments into the background.
    """

    start = time.time()
    # First, we convolve the background to reveal the points where we can
    # safely place a segment without overlapping them.
    valid = space.collisions(segment)
    convolution_duration = time.time() - start
    log(f"        (convolution took {convolution_duration:.3f} s)")

    if valid.size == 0:
        return None

    placements = []
    hits = 0
    tries = 0  # Number of times segment placement was unsuccesful.
    start = time.time()
    previously_selected = set()
    while hits < max_at_once:
        if tries >= max_tries:
            # Give up.
            log("  ! tries is exceeding max_tries")
            break
        if len(previously_selected) == len(valid[0]):
            log("  ! all open positions have been tried")
            break
        while True:
            selection = RNG.integers(0, len(valid[0]))
            if selection not in previously_selected:
                previously_selected.add(selection)
                break

        # Make sure that this placement does not overlap with another
        # previously selected location.
        location = placement_location(valid, selection, segment.voxels())
        prospect = np.where(segment.voxels()) + location[:, None]
        # Check for collisions at the prospective site.
        free = not np.any(
            space.squeezed_session_background[
                prospect[0, :], prospect[1, :], prospect[2, :]
            ]
        )

        if free:
            start = time.time()

            space.stamp(prospect)

            # If we are looking at a squeezed background, we need to correct the offset.
            # We also apply the segment's center translation.
            corrected_location = (
                location + space.squeezed_location_offset
            ) * space.resolution
            placements.append(corrected_location.tolist())
            hits += 1
        else:
            tries += 1
            log(
                f"        {tries = }/{max_tries},\thits = {hits}/{max_at_once}",
                end="\r",
            )
            placement_duration = time.time() - start
            if placement_duration * threshold_coefficient > convolution_duration:
                # The placement is taking half as long as a convolution
                # takes. At that point, it is probably cheaper to just run
                # another convolution.
                log(
                    f"  & placement_duration ({placement_duration:.6f}) * tc ({threshold_coefficient:.2f}) > convolution_duration"
                )
                break
    log("\x1b[K", end="")  # Clear the line since we used \r before.
    log(f"  . placed {hits} segments with {tries}/{max_tries} misses")

    return placements


def pack(
    space,
    segment,
    max_rotations,
    threshold_coefficient=1,
    max_at_once=10,
    max_tries=4,
) -> int:
    """
    Pack the background with a segment.
    """
    max_num = segment.target_number
    start_volume = np.sum(space.global_background)
    bgshape = space.global_background.shape
    background_volume = bgshape[0] * bgshape[1] * bgshape[2]  # TODO: np.product?

    segment_hits = 0
    rotations = 1  # The first round uses the identity rotation.

    # First, we try to just place the structure somewhere without the expensive convolution batch.
    while segment_hits < segment.target_number:
        remaining_segments = segment.target_number - segment_hits
        if remaining_segments > 100:
            target = remaining_segments // max_rotations
        else:
            target = remaining_segments

        start = time.time()
        placements = place(
            space.squeezed_global_background,
            space.squeezed_session_background,
            segment.voxels(),
            space.squeezed_location_offset,
            space.resolution,
            target,
            10000,  # Max tries for the naive placement.
            int(RNG.random() * 1000000),
        )
        duration = time.time() - start

        # FIXME: Make this a non-magic value.
        if len(placements) == 0 or len(placements) < target * 0.2:
            break

        n_placements = len(placements)
        segment_hits += n_placements
        segment.add_rotation(placements)
        segment.set_rotation(R.random(random_state=RNG).as_matrix())

        if n_placements != 0:
            log(
                f"%%%%%%% ({duration:.3f} s, {duration / n_placements:.6f} s per segment)"
            )
        else:
            log(f"%%%%%%% ({duration:.3f} s, no segments placed)")

    iteration = 0
    while segment_hits < segment.target_number:
        remaining_segments = max_num - segment_hits
        max_at_once_clamped = min(max_at_once, remaining_segments)
        assert max_at_once_clamped > 0
        log(
            f"  ^ trying to place {max_at_once_clamped} segments for this convolution",
        )
        print(
            f"    trying to place {max_at_once_clamped} segments for this convolution",
        )

        start = time.time()
        placements = place_batched(
            space,
            segment,
            max_at_once_clamped,
            max_tries,
            threshold_coefficient,
        )
        duration = time.time() - start

        if placements is None:
            log(
                f"    iteration {iteration} ended because the convolution produced no viable spots"
            )
            log(f"  * finished packing with {segment_hits} hits")
            break

        n_placements = len(placements)
        segment_hits += n_placements
        segment.add_rotation(placements)
        segment.set_rotation(R.random(random_state=RNG).as_matrix())

        if n_placements != 0:
            log(
                f"        ({duration:.3f} s, {duration / n_placements:.6f} s per segment)"
            )
        else:
            log(f"        ({duration:.3f} s, no segments placed)")

        if rotations >= max_rotations:
            log(f"        reached the maximum number of rotations ({max_rotations})")
            break

        rotations += 1
        iteration += 1

    end_volume = np.sum(space.global_background)
    print(
        f"  * finished packing with {segment_hits}/{max_num} hits ({start_volume / background_volume:.1%}->{end_volume / background_volume:.1%})"
    )
    return segment_hits


def format_placement_list(state):
    return json.dumps(
        {
            "title": state.title,
            "size": state.space.size,
            "topol_includes": state.topol_includes,
            "placements": [
                {
                    "name": segment.name,
                    "path": segment.path,
                    "batches": segment.batches,
                }
                for segment in state.segments
            ],
        }
    )


def setup_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="Pack a space.",
            prog="pack",
            epilog='"Maybe one minute is a bit optimistic."',
        )
    parser.add_argument(
        "path",
        type=pathlib.Path,
        help="The json configuration file to define the run.",
    )
    parser.add_argument(
        "--rearrange",
        action="store_true",
        help="Sort the input structures by approximate size to optimize packing.",
    )
    parser.add_argument(
        "--seed", default=None, type=int, help="Random number generator seed."
    )
    parser.add_argument(
        "--rotations",
        default=10,
        type=int,
        help="Set the number of random rotations per segment kind. (default: %(default)s)",
    )
    parser.add_argument(
        "--bead-radius",
        default=0.2,
        type=float,
        help="Set the bead radius that is considered during voxelization in nm. (default: %(default)s nm)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Use verbose output."
    )
    return parser


def configure(args=None):
    if args is None:
        parser = setup_parser()
        args = parser.parse_args()

    with open(args.path, "r") as config_file:
        config_src = config_file.read()

    config = Configuration(
        config_src,
        args.verbose,
        args.rearrange,
        args.seed,
        args.rotations,
        args.bead_radius,
    )

    global VERBOSE
    VERBOSE = config.verbose
    global RNG
    RNG = np.random.default_rng(seed=config.seed)

    return config


def main(state=None):
    if state is None:
        state = configure()
    space = state.space

    start = time.time()
    summary = []
    for i, segment in enumerate(state.segments):
        space.enter_session(segment.compartment_ids)
        segment_start = time.time()
        hits = pack(
            space,
            segment,
            state.rotations,
            # HACK: For now, this may be necessary.
            max_at_once=math.ceil(segment.target_number / state.rotations),
            max_tries=100,  # Maximum number of times to fail to place a segment.
        )
        segment_end = time.time()
        segment_duration = segment_end - segment_start
        space.exit_session()
        print(
            f"({i + 1}/{len(state.segments)}) Packed '{segment.name}' with a total of {hits} segments in {segment_duration:.3f} s.",
            flush=True,
        )
        nrots = len(segment.batches)
        summary.append(
            (segment.name, segment.target_number, hits, nrots, segment_duration)
        )
    end = time.time()
    packing_duration = end - start
    print(f"Packing process took {packing_duration:.3f} s.")

    placement_list_path = f"{state.output_dir}/{state.title}_placements.json"
    with open(placement_list_path, "w") as file:
        placement_list = format_placement_list(state)
        file.write(placement_list)
        print(f"Wrote placement list to '{placement_list_path}'.")

    print("Summary:")
    print("  idx \tname      \tnrots\ttarget\tplaced\ttime (s)\tremark")
    print("  ----\t----------\t-----\t------\t------\t--------\t------")
    nrots_tot = 0
    target_tot = 0
    hits_tot = 0
    for i, (name, target, hits, nrots, duration) in enumerate(summary):
        ok = " " if hits == target else "<"
        print(
            f"{i:>4}\t{name:>10}\t{nrots:>5}\t{target:>6}\t{hits:>6}\t{duration:8.2f}\t{ok}"
        )
        nrots_tot += nrots
        target_tot += target
        hits_tot += hits
    print(
        f"    \t          \t{nrots_tot:>5}\t{target_tot:>6}\t{hits_tot:>6}\t{packing_duration:8.2f}"
    )


if __name__ == "__main__":
    main()
