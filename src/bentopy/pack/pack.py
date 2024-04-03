import argparse
import json
import math
import time
import itertools
import pathlib

import numpy as np
from scipy.signal import fftconvolve
from scipy.spatial.transform import Rotation as R

from .config import Configuration

VERBOSE = False
ROTATIONS = 4
RNG = np.random.default_rng()


def log(*args, **kwargs):
    """
    Only print to console if verbose output is desired.
    """
    if VERBOSE:
        print(*args, **kwargs)


def trail_off(v: float) -> float:
    """
    Takes a value `v` in the range `0.0 <= v <= 1.0` and returns a probability
    according to a sigmoid density distribution, such that this probability
    is 1.0 for `v == 0.0` and goes to 0.0 as `v` approaches 1.0.
    """
    # See https://www.desmos.com/calculator/cdyqdxvjmy.
    return 1 - 1 / (1 + np.exp(-12 * (v - 0.5)))


def placement_location(valid, selection, segment):
    """
    Returns the valid voxel placement position for the segment.
    """
    segment_center = np.array(segment.shape) // 2
    valid_pos = valid[:, selection]
    return valid_pos - segment_center


def place(
    background,
    segment_voxels,
    max_at_once,
    max_tries,
    threshold_coefficient,
    resolution,
) -> list | None:
    """
    Place a number of segments into the background.
    """

    start = time.time()

    # FIXME(domains): Revisit and reformat this comment. Really needs to be correct ;)
    # We break up the background into a set of layers over the the longest axis.
    # We call these layers 'domains'. Each of these domains will be packed
    # sequentially, such that we place the segment over the entire space,
    # eventually. This breaking up into domains serves to reduce the peak
    # memory footprint during the convolution. For large spaces, the memory
    # required for a convolution of the entire space may exceed the available
    # amount. Hence, we want the ability to break up the space over which
    # the convolution is applied.
    # To make sure that the resulting packing remains well-stirred despite the
    # divided procedure, two passes are performed. The first we call the
    # 'white' pass, and it leaves clear division planes between the domains.
    # The second pass, which we call the 'black' pass, goes over the space with
    # an offset of half a domain size.
    # To ensure the placement densities are uniform overall, a placement
    # probability density function is used within a domain. It is shaped such
    # that the first half of the domain is a constant, full density. In this
    # first half, every valid placement that is considered will be accepted.
    # Within the second half of a domain, the placement probability density
    # will trail off to zero according to a sigmoid distribution. According to
    # this distribution, the closer a valid placement is to the end of the
    # domain, the less likely it becomes the placement is accepted.
    n_domains = 5
    domain_axis = np.argmax(background.shape)  # The longest axis.
    segments_per_domain = max_at_once / n_domains
    # The size of the background in the direction of the domain layers.
    background_size = background.shape[domain_axis]
    assert (
        background.shape[domain_axis] % n_domains == 0
    ), "For now, make sure we can neatly divide the background into domains"
    domain_size = background.shape[domain_axis] // n_domains
    assert (
        domain_size % 2 == 0
    ), "For now, make sure we can neatly divide a domain into two halves"
    half_domain_size = domain_size // 2
    # # We first pack the domains at the even indices, followed by the odd ones.
    # domain_indices = (*range(0, n_domains, 2), *range(1, n_domains, 2))
    white_domains = [
        ((2 * i) * half_domain_size, (2 * i + 2) * half_domain_size)
        for i in range(n_domains)
    ]
    # The black domains should also cover the first and last 'half' domains.
    black_domains = [
        (
            max(
                (2 * i - 1) * half_domain_size,
                0,
            ),
            min(
                (2 * i + 1) * half_domain_size,
                background_size,
            ),
        )
        for i in range(n_domains + 1)
    ]

    query = np.flip(segment_voxels)

    placements = []
    for domain_start, domain_end in itertools.chain(white_domains, black_domains):
        domain_slice = slice(domain_start, domain_end)
        if domain_axis == 0:
            domain = background[domain_slice, :, :]
        elif domain_axis == 1:
            domain = background[:, domain_slice, :]
        elif domain_axis == 2:
            domain = background[:, :, domain_slice]
        domain_offset = np.roll(np.array([domain_start, 0, 0]), domain_axis)
        print(f"{domain_offset = }")
        # First, we convolve the domain to reveal the points where we can
        # safely place a segment without overlapping them.
        collisions = fftconvolve(domain, query, mode="valid")
        valid_offset = np.array(query.shape) // 2
        convolution_duration = time.time() - start
        log(f"        (convolution took {convolution_duration:.3f} s)")

        # The valid placement points will have a value of 0. Since the floating
        # point operations leave some small errors laying around, we use a quite
        # generous cutoff.
        valid = np.array(np.where(collisions < 1e-4)) + valid_offset[:, None]
        if valid.size == 0:
            # If there are no valid placements, move on the next domain early.
            placements.append(None)
            continue

        domain_placements = []
        hits = 0
        tries = 0  # Number of times segment placement was unsuccesful.
        start = time.time()
        previously_selected = set()
        while hits < segments_per_domain and hits < max_at_once:
            if tries >= max_tries:
                # Give up.
                log("  ! tries is exceeding max_tries")
                break
            if len(previously_selected) == len(valid[0]):
                domain_placements = None
                break
            while True:
                selection = RNG.integers(0, len(valid[0]))
                if selection not in previously_selected:
                    previously_selected.add(selection)
                    break

            # Make sure that this placement does not overlap with another
            # previously selected location.
            location = placement_location(valid, selection, segment_voxels)
            prospect = np.where(segment_voxels) + location[:, None]
            # Check for collisions at the prospective site.
            free = not np.any(domain[prospect[0, :], prospect[1, :], prospect[2, :]])

            if free:
                start = time.time()

                domain_axis_location = location[domain_axis]
                # FIXME(domains): This is truly weird and can probably all be integrated in one well-formed if statement.
                is_half_domain = domain_end - domain_start == half_domain_size
                if is_half_domain:
                    # We need to deal with a half-sized domain at the start or end of the black pass.
                    if domain_start == 0:
                        acceptance_probability = trail_off(
                            domain_axis_location / half_domain_size
                        )
                    elif domain_end == background_size:
                        acceptance_probability = 1.0
                    else:
                        assert False, "Unreachable! Encountered a half domain that is not at the start or end of the background. This should be impossible."
                else:
                    # The common case.
                    if domain_axis_location < half_domain_size:
                        acceptance_probability = 1.0
                    if domain_axis_location >= half_domain_size:
                        acceptance_probability = trail_off(
                            (domain_axis_location - half_domain_size) / half_domain_size
                        )

                # Only actually place it if we accept it according to the
                # placement densitity distribution.
                if RNG.random() < acceptance_probability:
                    temp_selected_indices = prospect
                    domain[
                        temp_selected_indices[0, :],
                        temp_selected_indices[1, :],
                        temp_selected_indices[2, :],
                    ] = 1.0

                    domain_placements.append(
                        tuple(int(a) for a in (location + domain_offset) * resolution)
                    )
                    hits += 1
            else:
                tries += 1
                log(
                    f"        {tries = }/{max_tries},\thits = {hits}/{segments_per_domain}({max_at_once})",
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
        placements.append(domain_placements)

    print(f"placed {len(placements)} domains")

    # If the placements for each domain are all None, we need to let the caller
    # know. But if only some or none of them are None, we can just sum to
    # communicate the number of segments that were placed.
    if all(map(lambda domain_placements: domain_placements is None, placements)):
        print("None")
        return None
    else:
        r = []
        for domain_placements in placements:
            if domain_placements is not None:
                r.extend(domain_placements)
        print("return r")
        return r


def pack(
    background,
    segment,
    threshold_coefficient=1,
    max_iters=1000,
    max_at_once=10,
    max_tries=4,
    target_density=None,
) -> int:
    """
    Pack the background with a segment.
    """
    max_num = segment.target_number
    start_volume = np.sum(background)
    segment_volume = np.sum(segment)
    background_volume = (
        background.shape[0] * background.shape[1] * background.shape[2]
    )  # TODO: np.product?
    max_volume = background_volume - start_volume

    log("--> initiating packing")
    segment_hits = 0
    for iteration in range(max_iters):
        if segment_hits >= max_num:
            break
            log(f" -> iteration {iteration}")
        if target_density:
            # In case it could overshoot the target_density, limit the maximum
            # number of elements that can be placed per convolution.
            current_density = (np.sum(background) - start_volume) / max_volume
            remaining_density = target_density - current_density
            remaining_segment_volume = (
                remaining_density * max_volume
            )  # For this kind of segment.
            remaining_segments = math.ceil(remaining_segment_volume / segment_volume)
        else:
            remaining_segments = max_num - segment_hits
        max_at_once_clamped = min(max_at_once, remaining_segments)
        assert max_at_once_clamped > 0
        log(
            f"  ^ trying to place {max_at_once_clamped} segments for this convolution",
        )
        start = time.time()

        query = segment.voxels()
        placements = place(
            background,
            query,
            max_at_once_clamped,
            max_tries,
            threshold_coefficient,
            segment.resolution,
        )
        duration = time.time() - start
        log(f"        ({duration:.3f} s)")

        density = (np.sum(background) - start_volume) / max_volume
        if placements is None:
            log(
                f"    iteration {iteration} ended because the convolution produced no viable spots"
            )
            log(
                f"  * finished packing with {segment_hits} hits at a density of {density:.3f}"
            )
            break

        n_placements = len(placements)
        segment_hits += n_placements
        segment.add_rotation(placements)
        if n_placements != 0:
            log(f"        ({duration / n_placements:.6f} s per segment)")

        if target_density:
            # Do a density check.
            if density >= target_density:
                log(
                    f"    target density of {target_density:.2} was reached with {density:.2}"
                )
                break

        if target_density:
            log(
                f"    placed a total of {segment_hits}/{max_num} hits with a packing density of {density:.4f}/{target_density:.4f}"
            )
        else:
            log(
                f"    placed a total of {segment_hits}/{max_num} hits with a packing density of {density:.4f}"
            )

        segment.rotation = R.random(random_state=RNG).as_matrix()
    end_volume = np.sum(background)
    density = (end_volume - start_volume) / max_volume
    print(
        f"  * finished packing with {segment_hits}/{max_num} hits ({start_volume / background_volume:.1%}->{end_volume / background_volume:.1%}; ρ->{density:.3f})"
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
        "-v", "--verbose", action="store_true", help="Use verbose output."
    )
    return parser


def configure(args=None):
    if args is None:
        parser = setup_parser()
        args = parser.parse_args()

    with open(args.path, "r") as config_file:
        config_src = config_file.read()

    config = Configuration(config_src, args.verbose, args.rearrange, args.seed)

    global VERBOSE
    VERBOSE = config.verbose
    global RNG
    RNG = np.random.default_rng(seed=config.seed)

    return config


def main(state=None):
    if state is None:
        state = configure()
    background = state.space.background()

    start = time.time()
    for segment in state.segments:
        segment_background = state.space.background(
            segment.compartment_ids, onto=background
        )
        segment_start = time.time()
        hits = pack(
            segment_background,
            segment,
            # HACK: For now, this may be necessary.
            max_at_once=math.ceil(segment.target_number / 10),
            max_tries=100,  # Maximum number of times to fail to place a segment.
        )
        segment_end = time.time()
        segment_duration = segment_end - segment_start
        print(
            f"Packing '{segment.name}' with a total of {hits} segments took {segment_duration:.3f} s."
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
