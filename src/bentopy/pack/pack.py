import argparse
import json
import math
import os
import time
import pathlib

import MDAnalysis as MDA
import numpy as np
from scipy.signal import fftconvolve
from scipy.spatial.transform import Rotation as R

VERBOSE = False
ROTATIONS = 4
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
    padded_background,
    inside,  # A slice object to get the background from the padded_background.
    segment_voxels,
    max_at_once,
    max_tries,
    threshold_coefficient,
    resolution,
) -> list | None:
    """
    Place a number of segments into the background.
    """
    # First, we convolve the background to reveal the points where we can
    # safely place a segment without overlapping them.
    start = time.time()

    # The padded_background is padded very conservatively, and we don't need
    # all of its padding. We want to construct a slice object that we use to
    # select only the minimally needed section of the padded_background.
    minimal_padsize = np.array(segment_voxels.shape)
    # Minimal inside of the padded_background.
    inside_minimal = tuple(
        slice(sl.start - nps, sl.stop + nps)
        for (nps, sl) in zip(minimal_padsize, inside)
    )
    # The inside of collisions, such that its size equals the original background.
    inside_collisions = tuple(slice(p, -p) for p in minimal_padsize)

    query = np.flip(segment_voxels)
    # We pad the background in order to circumvent the edge effects of the
    # convolution. By padding and subsequently cropping the `collisions`
    # matrix, we make sure that there will be no out-of-bounds false positives
    # in the valid list.
    # TODO: There must be a more elegant way.
    collisions_padded = fftconvolve(
        padded_background[inside_minimal], query, mode="same"
    )
    convolution_duration = time.time() - start
    log(f"        (convolution took {convolution_duration:.3f} s)")

    # The valid placement points will have a value of 0. Since the floating
    # point operations leave some small errors laying around, we use a quite
    # generous cutoff.
    valid = np.array(np.where(collisions_padded[inside_collisions] < 1e-4))
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
            return None
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
        free = not np.any(
            padded_background[inside][prospect[0, :], prospect[1, :], prospect[2, :]]
        )

        if free:
            start = time.time()

            temp_selected_indices = prospect
            padded_background[inside][
                temp_selected_indices[0, :],
                temp_selected_indices[1, :],
                temp_selected_indices[2, :],
            ] = 1.0

            placements.append(tuple(int(a) for a in location * resolution))
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

    diagonal = segment.points().max(axis=1)
    padsize = int(np.ceil(np.linalg.norm(diagonal)))
    padded_background = np.pad(background, padsize, mode="constant", constant_values=2)
    # This slice object helps us retrieve the original background from a padded background array.
    inside = tuple(slice(padsize, -padsize) for _ in range(3))

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
            padded_background,
            inside,
            query,
            max_at_once_clamped,
            max_tries,
            threshold_coefficient,
            segment.resolution,
        )
        duration = time.time() - start
        log(f"        ({duration:.3f} s)")

        density = (np.sum(padded_background[inside]) - start_volume) / max_volume
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
    assert background.shape == padded_background[inside].shape
    background = padded_background[inside]
    end_volume = np.sum(background)
    density = (end_volume - start_volume) / max_volume
    print(
        f"  * finished packing with {segment_hits}/{max_num} hits ({start_volume / background_volume:.1%}->{end_volume / background_volume:.1%}; ρ->{density:.3f})"
    )
    return segment_hits


class Configuration:
    def __init__(self, json_src: str, verbose: bool, rearrange: bool, seed):
        config = json.loads(json_src)
        space = config["space"]
        self.space = Space(space["size"], space["resolution"], space["compartments"])
        segments = config["segments"]
        self.segments = [
            Segment(
                s["name"],
                s["number"],
                s["path"],
                self.space.resolution,
                s["compartments"],
            )
            for s in segments
        ]
        self.seed = seed
        if rearrange:
            if verbose:
                print("Order before:")
                for segment in self.segments:
                    print(f"\t{segment.name}")
            self.segments.sort(key=lambda seg: seg.voxels().sum(), reverse=True)
            if verbose:
                print("Order after:")
                for segment in self.segments:
                    print(f"\t{segment.name}")
            print("Rearranged the segments.")
        output = config["output"]
        self.title = output["title"]
        self.output_dir = output["dir"]
        if not os.path.isdir(self.output_dir):
            print(
                f"Output directory '{self.output_dir}' does not exist yet and will be created."
            )
            os.makedirs(self.output_dir)
        if "topol_includes" in output:
            self.topol_includes = output["topol_includes"]
        else:
            self.topol_includes = []
        self.verbose = verbose


class Voxels:
    def __init__(self, path):
        self.path = path
        self.voxels_cache = None
        self.mask_cache = None

    def voxels(self):
        if self.voxels_cache is None:
            if self.path.endswith(".npz"):
                with np.load(self.path) as npz:
                    self.voxels_cache = npz[npz.files[0]]  # FIXME: Seems silly.
            else:
                self.voxels_cache = np.load(self.path)
        return self.voxels_cache

    def mask(self, _width, _height, _depth, _padding):
        if self.mask_cache is None:
            voxels = self.voxels()
            # TODO: Apply the width, height, depth, padding parameters?
            self.mask_cache = np.where(voxels)
        return self.mask_cache


class Shape:
    kind_spherical = "spherical"
    kind_cuboid = "cuboid"
    kind_none = "none"
    available = [kind_spherical, kind_cuboid, kind_none]

    def __init__(self, kind):
        if kind not in self.available:
            raise ValueError(f"Unknown Shape kind: '{kind}' not in {self.available}")
        self.kind = kind
        self.mask_cache = None

    def mask(self, width, height, depth, padding):
        if self.mask_cache is None:
            if self.kind == self.kind_spherical:
                size = min(width, height, depth)
                # TODO: Make this elliptical or something idk.
                self.mask_cache = np.where(sphere_mask(size, size / 2 - padding))
            elif self.kind == self.kind_cuboid:
                self.mask_cache = cuboid_mask(width, height, depth, padding)
            elif self.kind == self.kind_none:
                self.mask_cache = tuple()
        return self.mask_cache


class Compartment:
    def __init__(self, id, shape=None, voxels=None):
        # Take care of the two incorrect input scenarios.
        both_none = shape is None and voxels is None
        both_some = shape is not None and voxels is not None
        if both_none or both_some:
            raise ValueError("Provide either a `shape` or `voxels` definition.")
        self.id = id
        if shape is not None:
            self.definition = Shape(shape)
        if voxels is not None:
            self.definition = Voxels(voxels["path"])

    def mask(self, width, height, depth, padding):
        return self.definition.mask(width, height, depth, padding)


class Space:
    def __init__(self, size, resolution, compartments):
        assert len(size) == 3, "The size of a space must be 3-dimensional"
        # Size in nm (unaffected by whatever value resolution has).
        self.size = size
        self.resolution = resolution
        self.compartments = []
        for c in compartments:
            id = c["id"]
            if "shape" in c:
                compartment = Compartment(id, shape=c["shape"])
            elif "voxels" in c:
                compartment = Compartment(id, voxels=c["voxels"])
            self.compartments.append(compartment)

    def background(self, compartment_ids=[], onto=None):
        # Adjust size and padding for resolution.
        size = (np.array(self.size) / self.resolution).astype(int)
        if onto is None:
            background = np.ones(size, dtype=np.float32)
        else:
            background = onto

        width, height, depth = size
        for compartment in filter(lambda c: c.id in compartment_ids, self.compartments):
            mask = compartment.mask(width, height, depth, padding=0)
            # Apply the mask
            background[mask] = 0
        return background


class Segment:
    def __init__(self, name, target_number, path, resolution, compartment_ids):
        self.name = name
        self.target_number = target_number
        self.path = path
        self.rotation = np.eye(3)  # We start out with the identity matrix.
        self.resolution = resolution
        self.compartment_ids = compartment_ids
        self._points = None
        # TODO: Perhaps this ought to be a dictionary, but then we would need a
        # method of normalizing rotation matrices to a lower precision in order
        # to make the bit patterns of ~identical rotations exactly the same.
        # Otherwise there'd be no benefit.
        self.batches = []  # A batch is a set of placements with a particular rotation.

    def points(self):
        """
        Returns the centered point cloud of the atoms in this Segment.

        Positions are given in nanometers.

        If no point cloud is stored internally, yet, it will be loaded from
        this Segment's path.
        """
        if self._points is None:
            self._points = load_points(self.path)

        return self._points

    def voxels(self):
        """
        Returns the voxels for this Segment's point cloud, with this
        Segment's rotation applied before voxelization.
        """
        return voxelize(self.rotation @ self.points().T, self.resolution, tighten=True)

    def add_rotation(self, placements):
        """
        Add a batch of placements with a particular rotation.

        The rotation is provided as an `[x, y, z]` list where the first value
        represents the rotation in degrees around the _x_-axis, etc.
        """
        self.batches.append((self.rotation.tolist(), placements))


def load_points(path):
    """
    Load atoms from a molecule file at `path` and return them as a point cloud.

    Positions in nanometers.
    """
    # TODO: Revisit the fact that we're creating the universe two times over
    # the runtime (depending on whether the placement list is rendered to a gro
    # file internally).
    # TODO: Consider whether it should be loaded into memory explicitly.
    u = MDA.Universe(path)
    # Convert from Å to nm.
    return u.atoms.positions / 10.0


def voxelize(points, resolution, tighten=False):
    """
    Return the voxels corresponding to a point cloud.

    The value of the provided resolution will be the final voxel size.
    """
    points = points / resolution
    mins = points.min(axis=1)
    points -= mins[:, None]

    maxs = np.ceil(points.max(axis=1)).astype(int)
    voxels = np.zeros((maxs[0], maxs[1], maxs[2]))  # FIXME: Better way?

    # TODO: This seems silly. Must be a better way.
    for x, y, z in points.T:
        voxels[int(x), int(y), int(z)] = 1

    if tighten:
        # And now, "tighten up those lines!"
        # TODO: We can do this in a more pretty manner, I'm sure.
        mins = np.min(np.where(voxels), axis=1)
        voxels = voxels[mins[0] :, mins[1] :, mins[2] :]
        maxs = np.max(np.where(voxels), axis=1) + 1
        voxels = voxels[: maxs[0], : maxs[1], : maxs[2]]

    return voxels


# TODO: There's got to be a nicer way of writing this function.
def sphere_mask(size, r):
    x = np.arange(0, size)
    y = np.arange(0, size)
    z = np.arange(0, size)
    cx = cy = cz = size / 2 - 1
    return (
        (x[:, None, None] - cx) ** 2
        + (y[None, :, None] - cy) ** 2
        + (z[None, None, :] - cz) ** 2
    ) < r**2


def cuboid_mask(width, height, depth, padding):
    np.s_[
        int(padding) : int(width) - int(padding),
        int(padding) : int(height) - int(padding),
        int(padding) : int(depth) - int(padding),
    ]


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
