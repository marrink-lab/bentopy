import json
import os

from .segment import Segment
from .space import Space

DEFAULT_COMPARTMENT_ID = "main"


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
        compartments = space.get("compartments")
        self.space = Space(
            space["size"],
            space["resolution"],
            compartments or [{"id": DEFAULT_COMPARTMENT_ID, "shape": "cuboid"}],
            space.get("constraint"),
        )
        segments = config["segments"]
        # Verify that if `output.compartments` is unset, the `compartments`
        # field is unset for every segment, as well.
        if compartments is None:
            assert all(s.get("compartments") is None for s in segments), """
            If no `compartments` are defined in the `space` section,
            compartment selections in segment definitions are not allowed.
            When no compartments are defined in the `space` section, a single
            default compartment is implied. To set each segment's compartment
            to the implied default compartment, remove the `compartments`
            field from all segment definitions."""
        self.segments = [
            Segment(
                s["name"],
                s["number"],
                s["path"],
                self.space.resolution,
                s.get("compartments") or [DEFAULT_COMPARTMENT_ID],
                s.get("rotation_axes"),
                s.get("center"),
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
