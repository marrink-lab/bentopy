import argparse
import json
import math
import time
import pathlib

import numpy as np
from scipy.spatial.transform import Rotation as R

from .config import Configuration

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


def place(
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
    max_iters=1000,
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
    for iteration in range(max_iters):
        if segment_hits >= max_num:
            break
        remaining_segments = max_num - segment_hits
        max_at_once_clamped = min(max_at_once, remaining_segments)
        assert max_at_once_clamped > 0
        log(
            f"  ^ trying to place {max_at_once_clamped} segments for this convolution",
        )
        start = time.time()

        placements = place(
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

        if n_placements != 0:
            log(f"        ({duration:.3f} s, {duration / n_placements:.6f} s per segment)")
        else:
            log(f"        ({duration:.3f} s, no segments placed)")

        if rotations >= max_rotations:
            log(f"        reached the maximum number of rotations ({max_rotations})")
            break

        segment.set_rotation(R.random(random_state=RNG).as_matrix())
        rotations += 1
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
        help="Set the number of random rotations per segment kind.",
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
        config_src, args.verbose, args.rearrange, args.seed, args.rotations
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
            f"({i + 1}/{len(state.segments)}) Packed '{segment.name}' with a total of {hits} segments in {segment_duration:.3f} s."
        )
    end = time.time()
    duration = end - start
    print(f"Packing process took {duration:.3f} s.")

    placement_list_path = f"{state.output_dir}/{state.title}_placements.json"
    with open(placement_list_path, "w") as file:
        placement_list = format_placement_list(state)
        file.write(placement_list)
        print(f"Wrote placement list to '{placement_list_path}'.")


if __name__ == "__main__":
    main()
