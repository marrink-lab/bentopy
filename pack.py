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


def placement_location(valid, selection, segment):
    """
    Returns the valid voxel placement position for the segment.
    """
    segment_shape = np.array(segment.shape)
    return (
        np.array([valid[0][selection], valid[1][selection]])
        - (segment_shape // 2 + segment_shape % 2)
        + 1
    )


def place_segment_convolve(
    background,
    segment_voxels,
    max_at_once,
    max_tries,
    threshold_coefficient,
    resolution,
) -> list | None:
    query = np.flip(segment_voxels, axis=(0, 1))
    start = time.time()
    collisions = fftconvolve(background, query, mode="same")
    convolution_duration = time.time() - start
    if VERBOSE:
        print(f"        (convolution took {convolution_duration:.3f} s)")

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
            if VERBOSE:
                print("  ! tries is exceeding max_tries")
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
        prospect = (
            np.where(segment_voxels) + location[:, None] - np.array([0, 1])[:, None]
        )
        # Check for collisions at the prospective site.
        free = not np.any(background[*prospect])

        if free:
            start = time.time()

            temp_selected_indices = prospect
            background[*temp_selected_indices] = 1.0

            placements.append((*(int(a) for a in location * resolution), 0.0))
            hits += 1
        else:
            tries += 1
            if VERBOSE and tries % 10 == 0:
                print(
                    f"        {tries = }/{max_tries},\thits = {hits}/{max_at_once}",
                    end="\r",
                )
            placement_duration = time.time() - start
            if placement_duration * threshold_coefficient > convolution_duration:
                # The placement is taking half as long as a convolution
                # takes. At that point, it is probably cheaper to just run
                # another convolution.
                if VERBOSE:
                    print(
                        f"  & placement_duration ({placement_duration:.6f}) * tc ({threshold_coefficient:.2f}) > convolution_duration"
                    )
                break
    if VERBOSE:
        print("\x1b[K", end="")  # We used a \r before.
    if VERBOSE:
        print(f"  . placed {hits} segments with {tries}/{max_tries} misses")

    return placements


def fill_background(
    background,
    segment,
    resolution,
    threshold_coefficient=1,
    max_iters=1000,
    max_at_once=10,
    max_tries=4,
    target_density=None,
) -> int:
    max_num = segment.target_number
    start_volume = np.sum(background)
    segment_volume = np.sum(segment)
    background_volume = background.shape[0] * background.shape[1]  # TODO: np.product?
    max_volume = background_volume - start_volume
    if VERBOSE:
        print("--> initiating packing")
    rotation = [0, 0, 0]
    segment_hits = 0
    for iteration in range(max_iters):
        if segment_hits >= max_num:
            break
        if VERBOSE:
            print(f" -> iteration {iteration}")
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
        if VERBOSE:
            print(
                f"  ^ trying to place {max_at_once_clamped} segments for this convolution",
            )
        start = time.time()
        placements = place_segment_convolve(
            background,
            segment.voxels(),
            max_at_once_clamped,
            max_tries,
            threshold_coefficient,
            resolution,
        )
        duration = time.time() - start
        if VERBOSE:
            print(f"        ({duration:.3f} s)")

        if placements is None:
            if VERBOSE:
                print(
                    f"    iteration {iteration} ended because the convolution produced no viable spots"
                )
            density = (np.sum(background) - start_volume) / max_volume
            if VERBOSE:
                print(
                    f"  * finished packing with {segment_hits} hits at a density of {density:.3f}"
                )
            break

        n_placements = len(placements)
        segment_hits += n_placements
        segment.add_rotation(rotation, placements)
        if VERBOSE and n_placements != 0:
            print(f"        ({duration / n_placements:.6f} s per segment)")

        density = (np.sum(background) - start_volume) / max_volume
        if target_density:
            # Do a density check.
            if density >= target_density:
                if VERBOSE:
                    print(
                        f"    target density of {target_density:.2} was reached with {density:.2}"
                    )
                break

        if target_density:
            if VERBOSE:
                print(
                    f"    placed a total of {segment_hits}/{max_num} hits with a packing density of {density:.4f}/{target_density:.4f}"
                )
        else:
            if VERBOSE:
                print(
                    f"    placed a total of {segment_hits}/{max_num} hits with a packing density of {density:.4f}"
                )

        # Silly rotation for the next round while we're testing.
        segment._voxels = np.rot90(segment.voxels())
        rotation[2] += 90  # Update the z-axis to reflect the rotation.
    end_volume = np.sum(background)
    density = (end_volume - start_volume) / max_volume
    print(
        f"  * finished packing with {segment_hits}/{max_num} hits ({start_volume / background_volume:.1%}->{end_volume / background_volume:.1%}; ρ->{density:.3f})"
    )
    return segment_hits


class Configuration:
    def __init__(self, json_src: str, verbose: bool):
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
        self.size = size[:2]  # For now, we do a lil 2D thing.
        self.resolution = resolution
        self.shape = shape
        # Padding in nm (unaffected by whatever value resolution has).
        self.padding = padding

    def background(self):
        # Adjust size and padding for resolution.
        size = (np.array(self.size) / self.resolution).astype(int)
        padding = self.padding / self.resolution

        background = np.ones(size, dtype=np.float32)
        width, height = size
        if VERBOSE:
            print(f"{self.shape = }")
        if self.shape == "circular":
            size = min(width, height)
            # TODO: Make this elliptical or something idk.
            mask = circle_mask(size, size // 2 - padding)
        elif self.shape == "rectangular":
            mask = np.index_exp[
                padding : width - padding,
                padding : height - padding,
            ]
        else:
            raise ValueError
        background[mask] = 0
        return background


class Segment:
    def __init__(self, name, target_number, path, resolution):
        self.name = name
        self.target_number = target_number
        self.path = path
        self._resolution = resolution
        self._voxels = None
        self.batches = []  # A batch is a set of placements with a particular rotation.

    def voxels(self):
        if self._voxels is None:
            # FIXME: We just take one slice of the voxels cloud, for now. Revisit for 3D.
            voxels = structure_to_3d(self.path, self._resolution)
            zmid = voxels.shape[2] // 2
            # And now, "tighten up those lines!"
            slice = voxels[:, :, zmid]
            mins = np.min(np.where(slice), axis=1)
            slice = slice[mins[0] :, mins[1] :]
            maxs = np.max(np.where(slice), axis=1) + 1
            slice = slice[: maxs[0], : maxs[1]]
            self._voxels = slice
        return self._voxels

    def add_rotation(self, rotation, placements):
        """
        Add a batch of placements with a particular rotation.

        The rotation is provided as an `[x, y, z]` list where the first value
        represents the rotation in degrees around the _x_-axis, etc.
        """
        # Convert the rotation to a 3×3 rotation matrix.
        r = R.from_euler("zyx", rotation, degrees=True).as_matrix().tolist()
        self.batches.append((r, placements))


def plot_voxels(voxels):
    ax = plt.figure().add_subplot(projection="3d")
    ax.voxels(voxels, linewidth=0.5)
    ax.set(xlabel="r", ylabel="g", zlabel="b")
    ax.set_aspect("equal")
    plt.show()


def circle_mask(size, r):
    x = np.arange(0, size)
    y = np.arange(0, size)
    cx = size // 2
    cy = size // 2
    return (x[None, :] - cx) ** 2 + (y[:, None] - cy) ** 2 < r**2


def structure_to_3d(path, resolution):
    # TODO: Revisit the fact that we're creating the universe two times over the runtime (depending on the output format).
    u = MDA.Universe(path)
    # Convert from Å to nm and adjust for resulotion.
    positions = u.atoms.positions / 10.0 / resolution
    mins = positions.min(axis=0)
    positions -= mins

    maxs = np.ceil(positions.max(axis=0)).astype(int)
    voxels = np.zeros((maxs[0], maxs[1], maxs[2]))

    for x, y, z in positions:
        voxels[int(x), int(y), int(z)] = 1

    # plot_voxels(voxels)
    return voxels


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
        "-v", "--verbose", action="store_true", help="Use verbose output"
    )
    args = parser.parse_args()

    config_path = args.config[0]
    with open(config_path, "r") as config_file:
        config_src = config_file.read()
    return Configuration(config_src, args.verbose)


def render_to_gro(path, segments, box):
    start_file_writing = time.time()
    with open(path, "w") as gro:
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
                mins = np.min(rotated_positions, axis=0)
                # FIXME: The indexing here is for 2D, to get nice pancakes.
                rotated_positions[:, :2] -= mins[:2]
                for idx, placement in enumerate(placements):
                    dx, dy, dz = placement
                    # This 0.0 will become dz when we go 3D.
                    translation = np.array((dx, dy, 0.0))
                    positions = rotated_positions + translation
                    # TODO: Select relevant parts only.
                    for prefix, position in zip(prefixes, positions):
                        posx, posy, posz = position
                        line = f"{prefix}{posx:8.3f}{posy:8.3f}{posz:8.3f}"
                        gro.write(line + "\n")
                        total_atoms += 1

        # v1x, v2y, v3z = np.max(positions, axis=1)
        if len(box) == 2:
            v1x, v2y, v3z = *box, 10.0  # HACK: This is rather temporary.
        elif len(box) == 3:
            v1x, v2y, v3z = box
        box_vectors = f"{v1x:.3f} {v2y:.3f} {v3z:.3f}"
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
            config.space.resolution,
            max_at_once=math.ceil(
                segment.target_number / ROTATIONS  # Because we do want rotations!
            ),  # ???: Figure out where to go from here. May become a 'pretty parameter'.
            max_tries=400,  # Maximum number of times to fail to place a segment.
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
                    "size": [*config.space.size, 0],  # HACK: This will become 3D later.
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
        path = f"{config.output_dir}/{config.title}.png"
        plt.imsave(path, background)
        print(f"Wrote debug image to '{path}'.")


if __name__ == "__main__":
    main()
