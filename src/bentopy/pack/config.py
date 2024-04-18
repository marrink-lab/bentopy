import json
import os

from .segment import Segment
from .space import Space


class Configuration:
    def __init__(
        self,
        json_src: str,
        verbose: bool,
        rearrange: bool,
        seed: int,
        rotations: int,
    ):
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
        self.rotations = rotations
