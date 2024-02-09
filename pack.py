import argparse
import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import fftconvolve
import MDAnalysis as MDA
from MDAnalysis import transformations

VERBOSE = False
ROTATIONS = 4


def bounding_box(arr):
    rows = np.any(arr, axis=1)
    cols = np.any(arr, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax


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
    background, segment, max_at_once, max_tries, threshold_coefficient
) -> int:
    query = np.flip(segment, axis=(0, 1))
    start = time.time()
    collisions = fftconvolve(background, query, mode="same")
    # plt.imshow(collisions)
    # plt.show()
    convolution_duration = time.time() - start
    if VERBOSE:
        print(f"        (convolution took {convolution_duration:.3f} s)")

    valid = np.where(collisions < 1e-4)
    if valid[0].size > 0:
        single = False
        if single:
            # selection = 0 # Important for consistent benchmarking.
            selection = np.random.randint(0, len(valid[0]))

            # For now, just place once per convolution.
            locations = [placement_location(valid, selection, segment)]
        else:
            # The worst-case segment diameter.
            segment_diameter = np.linalg.norm(np.array(segment.shape))
            locations = []
            tries = 0  # Number of times segment placement was unsuccesful.
            start = time.time()
            while len(locations) < max_at_once:
                if tries >= max_tries:
                    # Give up.
                    if VERBOSE:
                        print("  ! tries is exceeding max_tries")
                    break
                selection = np.random.randint(0, len(valid[0]))

                # Check whether this selection is not too close to another
                # previously selected location.
                free = True
                for location in locations:
                    attempt = placement_location(valid, selection, segment)
                    distance = np.linalg.norm(attempt - location)
                    if distance < segment_diameter:
                        free = False
                        break

                if free:
                    start = time.time()
                    locations.append(placement_location(valid, selection, segment))
                else:
                    tries += 1
                    if tries % 10 == 0:
                        if VERBOSE:
                            print(
                                f"        {tries = }/{max_tries},\thits = {len(locations)}/{max_at_once}",
                                end="\r",
                            )
                    placement_duration = time.time() - start
                    if (
                        placement_duration * threshold_coefficient
                        > convolution_duration
                    ):
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
                print(
                    f"  . placed {len(locations)} segments with {tries}/{max_tries} misses"
                )

        # ??? May want to move this elsewhere.
        selected_voxel_indices = np.asarray(np.where(segment))

        hits = []
        for location in locations:
            temp_selected_indices = selected_voxel_indices.copy()
            for dim in range(len(selected_voxel_indices)):
                temp_selected_indices[dim] += location[dim]

            background[*temp_selected_indices] = 1
            hits.append(tuple(location))

        return hits
    return []


def fill_background(
    background,
    segment,
    max_num,
    threshold_coefficient=25,
    max_iters=10000,
    max_at_once=10,
    max_tries=4,
    target_density=None,
) -> list:
    start_volume = np.sum(background)
    segment_volume = np.sum(segment)
    background_volume = background.shape[0] * background.shape[1]  # TODO: np.product?
    max_volume = background_volume - start_volume
    if VERBOSE:
        print("--> initiating packing")
    segment = np.copy(segment)
    rotation_state = 0
    segment_hits = 0
    segment_placements = [
        [] for _ in range(ROTATIONS)
    ]  # A list of placements corresponding to each rotation.
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
        hits = place_segment_convolve(
            background, segment, max_at_once_clamped, max_tries, threshold_coefficient
        )
        duration = time.time() - start
        if VERBOSE:
            print(f"        ({duration:.3f} s)")
        segment_hits += len(hits)

        if len(hits) == 0:
            if VERBOSE:
                print(f"    iteration {iteration} ended with 0 hits")
            density = (np.sum(background) - start_volume) / max_volume
            if VERBOSE:
                print(
                    f"  * finished packing with {segment_hits} hits at a density of {density:.3f}"
                )
            break

        if VERBOSE:
            print(f"        ({duration / hits:.6f} s per segment)")
        segment_placements[rotation_state].extend(hits)

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
        segment = np.rot90(segment)
        rotation_state = (rotation_state + 1) % ROTATIONS
    end_volume = np.sum(background)
    density = (end_volume - start_volume) / max_volume
    print(
        f"  * finished packing with {segment_hits}/{max_num} hits ({start_volume / background_volume:.1%}->{end_volume / background_volume:.1%}; ρ->{density:.3f})"
    )
    return segment_placements


class Configuration:
    def __init__(self, json_src: str, verbose: bool):
        config = json.loads(json_src)
        space = config["space"]
        mask = space["mask"]
        self.space = Space(space["size"], mask["shape"], mask["padding"])
        segments = config["segments"]
        self.segments = [Segment(s["name"], s["number"], s["path"]) for s in segments]
        output = config["output"]
        self.output_path = output["path"]
        self.output_config = output
        self.verbose = verbose


class Space:
    def __init__(self, size, shape, padding):
        assert len(size) == 3, "The size of a space must be 3-dimensional"
        self.size = size[:2] # For now, we do a lil 2D thing.
        self.shape = shape
        self.padding = padding

    def background(self):
        # The playing field.
        background = np.ones(self.size, dtype=np.float32)
        width, height = self.size
        if VERBOSE:
            print(f"{self.shape = }")
        if self.shape == "circular":
            size = min(width, height)
            # TODO: Make this elliptical or something idk.
            mask = circle_mask(size, size // 2 - self.padding)
        elif self.shape == "rectangular":
            mask = np.index_exp[
                self.padding : width - self.padding,
                self.padding : height - self.padding,
            ]
        else:
            raise ValueError
        background[mask] = 0
        return background


class Segment:
    def __init__(self, name, number, path):
        self.name = name
        self.number = number
        self.path = path
        self._voxels = None

    def voxels(self):
        if not self._voxels:
            self._voxels = structure_to_2d(self.path)
        return self._voxels


def circle_mask(size, r):
    x = np.arange(0, size)
    y = np.arange(0, size)
    cx = size // 2
    cy = size // 2
    return (x[None, :] - cx) ** 2 + (y[:, None] - cy) ** 2 < r**2


def structure_to_2d(path):
    # TODO: Revisit the fact that we're creating the universe two times over the runtime (depending on the output format).
    u = MDA.Universe(path)
    positions = u.atoms.positions / 10.0  # Convert from Å to nm.
    center = positions.mean(axis=0)
    positions -= center
    slicer = positions[:, 2] < 0.4
    slice = positions[slicer][:, :2]

    xmin = slice[:, 0].min()
    ymin = slice[:, 1].min()
    slice[:, 0] -= xmin
    slice[:, 1] -= ymin
    xmax = slice[:, 0].max()
    ymax = slice[:, 1].max()
    arr = np.zeros((int(xmax) + 1, int(ymax) + 1))

    for x, y in slice:
        arr[int(x), int(y)] = 1

    return arr


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
        maxx, maxy, maxz = 0.0, 0.0, 0.0
        title = "pack"
        print(title, file=gro)  # Title.
        print("------------", file=gro)  # Placeholder for the number of atoms.
        total_atoms = 0
        for resnum, segment in enumerate(segments):
            u = MDA.Universe(segment.path)

            prefixes = []
            for atom in u.atoms:
                # resnum = atom.resnum
                resname = atom.resname
                atomname = atom.name
                atomnum = atom.ix
                # FIXME: resnum should not not actually be the same for every one of them I think.
                # This sequence is constant between
                prefix = f"{resnum:5d}{resname:<5s}{atomname:>5s}{atomnum:5d}"
                prefixes.append(prefix)

            ts = u.trajectory.ts
            ts.positions /= 10.0  # Convert from Å to nm.
            for rotation, rotation_group in enumerate(segment.placements):
                # TODO: This ts bs can probably go at some point. But for now we're going for a proof-of-concept.
                angle = 90 * rotation
                ag = u.atoms
                d = [0, 0, 1]  # Rotate around the z axis.
                rotated_positions = transformations.rotate.rotateby(
                    angle, direction=d, ag=ag
                )(ts).positions

                center = rotated_positions.mean(axis=0)
                rotated_positions -= center
                for idx, placement in enumerate(rotation_group):
                    dx, dy = placement
                    translation = np.array((dx, dy, 0.0))
                    positions = rotated_positions + translation
                    # TODO: Select relevant parts only.
                    for prefix, position in zip(prefixes, positions):
                        posx, posy, posz = position
                        line = f"{prefix}{posx:8.3f}{posy:8.3f}{posz:8.3f}"
                        print(line, file=gro)
                        total_atoms += 1

        # v1x, v2y, v3z = np.max(positions, axis=1)
        if len(box) == 2:
            v1x, v2y, v3z = *box, 10.0 # HACK: This is rather temporary.
        elif len(box) ==3:
            v1x, v2y, v3z = box
        box_vectors = f"{v1x:.6f} {v2y:.6f} {v3z:.6f}"
        print(box_vectors, file=gro)

        # Go back to the start and write the number of atoms.
        gro.seek(len(title) + 1)
        print(f"{total_atoms:>12d}", file=gro)


def main():
    config = configure()
    background = config.space.background()
    global VERBOSE
    VERBOSE = config.verbose

    start = time.time()
    for segment in config.segments:
        segment_start = time.time()
        segment_placements = fill_background(
            background,
            segment.voxels(),
            segment.number,
            max_at_once=math.ceil(
                segment.number / ROTATIONS  # Because we do want rotations!
            ),  # ???: Figure out where to go from here. May become a 'pretty parameter'.
            max_tries=4,  # Maximum number of times to fail to place a segment.
        )
        segment.placements = segment_placements
        segment_end = time.time()
        segment_duration = segment_end - segment_start
        print(f"packing '{segment.name}' took {segment_duration:.3f} s")
    end = time.time()
    duration = end - start
    print(f"packing process took {duration:.3f} s")

    # TODO: This is incorrect since the output config is different than this is.
    if config.output_path.endswith(".gro"):
        render_to_gro(config.output_path, config.segments, config.space.size)
    elif config.output_path.endswith(".pdb"):
        render_to_pdb(config.output_path, config.segments, config.space.size)

    # TODO: Remove when there's proper output or at least make it depend on debug output.
    plt.imsave(f"{config.output_path}.png", background)
    plt.show()


if __name__ == "__main__":
    main()
