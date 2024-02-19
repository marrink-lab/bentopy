import argparse
import json
import math
import os
import time

import matplotlib.pyplot as plt
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
    valid_pos = np.array(
        (valid[0][selection], valid[1][selection], valid[2][selection])
    )
    return valid_pos - segment_center


def place_segment_convolve(
    background,
    segment_voxels,
    max_at_once,
    max_tries,
    threshold_coefficient,
    resolution,
) -> list | None:
    # First, we convolve the background to reveal the points where we can
    # safely place a segment without overlapping them.
    start = time.time()
    query = np.flip(segment_voxels)
    # We pad the background in order to circumvent the edge effects of the
    # convolution. By padding and subsequently cropping the `collisions`
    # matrix, we make sure that there will be no out-of-bounds false positives
    # in the valid list.
    # TODO: The choice of pad_width could be optimized, but I doubt it would
    # make any significant difference.
    padwidth = max(np.array(query.shape))
    padded_background = np.pad(background, padwidth, mode="constant", constant_values=2)
    # TODO: There must be a more elegant way.
    collisions = fftconvolve(padded_background, query, mode="same")[
        padwidth:-padwidth, padwidth:-padwidth, padwidth:-padwidth
    ]
    # TODO: Maybe we can just remove this later.
    assert collisions.shape == background.shape
    convolution_duration = time.time() - start
    log(f"        (convolution took {convolution_duration:.3f} s)")

    # The valid placement points will have a value of 0. Since the floating
    # point operations leave some small errors laying around, we use a quite
    # generous cutoff.
    valid = np.where(collisions < 1e-4)
    valid_spots = valid[0].size
    if valid_spots == 0:
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
            selection = np.random.randint(0, len(valid[0]))
            if selection not in previously_selected:
                previously_selected.add(selection)
                break

        # Make sure that this placement does not overlap with another
        # previously selected location.
        location = placement_location(valid, selection, segment_voxels)
        prospect = np.where(segment_voxels) + location[:, None]
        # Check for collisions at the prospective site.
        free = not np.any(background[*prospect])

        if free:
            start = time.time()

            temp_selected_indices = prospect
            background[*temp_selected_indices] = 1.0

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


def fill_background(
    background,
    segment,
    threshold_coefficient=1,
    max_iters=1000,
    max_at_once=10,
    max_tries=4,
    target_density=None,
) -> int:
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
        placements = place_segment_convolve(
            background,
            segment.voxels(),
            max_at_once_clamped,
            max_tries,
            threshold_coefficient,
            segment.resolution,
        )
        duration = time.time() - start
        log(f"        ({duration:.3f} s)")

        if placements is None:
            log(
                f"    iteration {iteration} ended because the convolution produced no viable spots"
            )
            density = (np.sum(background) - start_volume) / max_volume
            log(
                f"  * finished packing with {segment_hits} hits at a density of {density:.3f}"
            )
            break

        n_placements = len(placements)
        segment_hits += n_placements
        segment.add_rotation(placements)
        if n_placements != 0:
            log(f"        ({duration / n_placements:.6f} s per segment)")

        density = (np.sum(background) - start_volume) / max_volume
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


class Configuration:
    def __init__(self, json_src: str, verbose: bool, rearrange: bool):
        config = json.loads(json_src)
        space = config["space"]
        mask = space["mask"]
        self.space = Space(
            space["size"], space["resolution"], mask["shape"], mask["padding"]
        )
        segments = config["segments"]
        self.segments = [
            Segment(s["name"], s["number"], s["path"], self.space.resolution)
            for s in segments
        ]
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
            os.makedir(self.output_dir)
        self.render = output["render"]
        if "debug_image" in output:
            self.debug_image = output["debug_image"]
        else:
            self.debug_image = None
        self.output_placement_list = output["placement_list"]
        self.verbose = verbose


class Space:
    def __init__(self, size, resolution, shape, padding):
        assert len(size) == 3, "The size of a space must be 3-dimensional"
        # Size in nm (unaffected by whatever value resolution has).
        self.size = size
        self.resolution = resolution
        self.shape = shape
        # Padding in nm (unaffected by whatever value resolution has).
        self.padding = padding

    def background(self):
        # Adjust size and padding for resolution.
        size = (np.array(self.size) / self.resolution).astype(int)
        padding = self.padding / self.resolution

        background = np.ones(size, dtype=np.float32)
        width, height, depth = size
        log(f"{self.shape = }")
        if self.shape == "circular":
            size = min(width, height, depth)
            # TODO: Make this elliptical or something idk.
            mask = sphere_mask(size, size // 2 - padding)
        elif self.shape == "rectangular":
            mask = np.index_exp[
                padding : width - padding,
                padding : height - padding,
                padding : depth - padding,
            ]
        elif self.shape == "none":
            mask = tuple()
        else:
            raise ValueError
        background[mask] = 0
        return background


class Segment:
    def __init__(self, name, target_number, path, resolution):
        self.name = name
        self.target_number = target_number
        self.path = path
        self.rotation = np.eye(3)  # We start out with the identity matrix.
        self.resolution = resolution
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


def plot_voxels(voxels):
    ax = plt.figure().add_subplot(projection="3d")
    ax.voxels(voxels, linewidth=0.5)
    ax.set(xlabel="r", ylabel="g", zlabel="b")
    ax.set_aspect("equal")
    plt.show()


# TODO: There's got to be a nicer way of writing this function.
def sphere_mask(size, r):
    x = np.arange(0, size)
    y = np.arange(0, size)
    z = np.arange(0, size)
    cx = size // 2
    cy = size // 2
    cz = size // 2
    return
    (
        (x[None, :, :] - cx) ** 2
        + (y[:, None, :] - cy) ** 2
        + (z[:, :, None] - cz) ** 2
    ) < r**2


def configure():
    parser = argparse.ArgumentParser(
        description="Pack a space.",
        prog="pack",
        epilog='"Maybe one minute is a bit optimistic."',
    )
    parser.add_argument(
        "config", metavar="path", nargs=1, type=str, help="json file to define the run"
    )
    parser.add_argument(
        "--rearrange",
        action="store_true",
        help="Sort the input structures by approximate size to optimize packing",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Use verbose output"
    )
    args = parser.parse_args()

    config_path = args.config[0]
    with open(config_path, "r") as config_file:
        config_src = config_file.read()
    return Configuration(config_src, args.verbose, args.rearrange)


def render_to_gro(path, segments, box):
    start_file_writing = time.time()
    with open(path, "w") as gro:
        # Format box vecs now, to prevent discovering if they're bad much later.
        v1x, v2y, v3z = box
        box_vectors = f"{v1x:.3f} {v2y:.3f} {v3z:.3f}"

        title = "pack"
        gro.write(title + "\n")
        gro.write("------------" + "\n")  # Placeholder for the number of atoms.
        total_atoms = 0
        for resnum, segment in enumerate(segments):
            u = MDA.Universe(segment.path)

            prefixes = []
            for atom in u.atoms:
                # resnum = atom.resnum
                resname = atom.resname
                atomname = atom.name
                atomnum = atom.ix + 1
                # FIXME: resnum should not not actually be the same for every one of them I think.
                # This sequence is constant between
                prefix = f"{resnum:>5d}{resname:<5s}{atomname:>5s}{atomnum:5d}"
                prefixes.append(prefix)

            ts = u.trajectory.ts
            ts.positions /= 10.0  # Convert from Å to nm.
            center = ts.positions.mean(axis=0)
            ts.positions -= center
            for rotation, placements in segment.batches:
                r = R.from_matrix(rotation)
                rotated_positions = r.apply(ts.positions)
                for idx, placement in enumerate(placements):
                    dx, dy, dz = placement
                    translation = np.array((dx, dy, dz))
                    positions = rotated_positions + translation
                    # TODO: Select relevant parts only.
                    for prefix, position in zip(prefixes, positions):
                        posx, posy, posz = position
                        line = f"{prefix}{posx:8.3f}{posy:8.3f}{posz:8.3f}"
                        gro.write(line + "\n")
                        total_atoms += 1

        gro.write(box_vectors + "\n")

        # Go back to the start and write the number of atoms.
        gro.seek(len(title) + 1)
        gro.write(f"{total_atoms:>12d}\n")
    end_file_writing = time.time()
    print(f"Wrote '{path}' in {end_file_writing - start_file_writing:.3f} s.")


def main():
    config = configure()
    background = config.space.background()
    global VERBOSE
    VERBOSE = config.verbose

    start = time.time()
    for segment in config.segments:
        segment_start = time.time()
        hits = fill_background(
            background,
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

    if config.output_placement_list:
        placement_list_path = f"{config.output_dir}/{config.title}_placements.json"
        with open(placement_list_path, "w") as outfile:
            placement_list_dump = json.dumps(
                {
                    "size": config.space.size,
                    "placements": [
                        {
                            "name": segment.name,
                            "path": segment.path,
                            "batches": segment.batches,
                        }
                        for segment in config.segments
                    ],
                }
            )
            outfile.write(placement_list_dump)
            print(f"Wrote placement list to '{placement_list_path}'.")

    # TODO: This is incorrect since the output config is different than this is.
    if config.render:
        render_to_gro(
            f"{config.output_dir}/{config.title}.gro",
            config.segments,
            config.space.size,
        )

    if config.debug_image:
        _path = f"{config.output_dir}/{config.title}.png"
        # plt.imsave(path, background)
        # print(f"Wrote debug image to '{path}'.")
        print("ERROR: CANNOT WRITE IMAGE TO PATH")


if __name__ == "__main__":
    main()
