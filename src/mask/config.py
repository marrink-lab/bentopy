import argparse
from pathlib import Path

EPILOG = """
Compartments are determined using the MDVContainment package created by Bart Bruininks.
More information is available at <https://github.com/BartBruininks/mdvcontainment>.
"""
DESCRIPTION = "Set up masks based on structures and compartment segmentations."

DEFAULT_CONTAINMENT_RESOLUTION = 1.0
DEFAULT_MASK_RESOLUTION = 0.5


def process_labels(labels: str) -> list[int]:
    return [int(label) for label in labels.split(",")]


def setup_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
            prog="mask",
            epilog=EPILOG,
        )
    else:
        parser.description = DESCRIPTION
        parser.epilog = EPILOG
    parser.add_argument(
        "input",
        type=Path,
        help="Input structure file to subject to segmentation.",
    )
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="Output path for the resulting voxel mask.",
    )
    parser.add_argument(
        "--no-interactive",
        dest="interactive",
        action="store_false",
        help="""Do not ask for command line input at runtime. May be desirable 
        when all inputs are known and in a scripted context.
        (normally interactive)""",
    )
    parser.add_argument(
        "--closing",
        action="store_true",
        help="""Use binary closing to fill small holes in compartments. For analysis
        of CG structures with a voxel resolution <1.0 nm this is highly recommended.
        """,
    )
    parser.add_argument(
        "--slab",
        action="store_true",
        help="""Determine the containment according to slab periodicity.
        This setting is suitable for simple systems where an 'outside' region can be trivially determined.
        When enabled, it may speed up containment determination.
        """,
    )
    parser.add_argument(
        "--containment-resolution",
        default=DEFAULT_CONTAINMENT_RESOLUTION,
        type=float,
        help="Resolution for the compartment finding routine (nm). (default: %(default)s nm)",
    )
    parser.add_argument(
        "--mask-resolution",
        default=DEFAULT_MASK_RESOLUTION,
        type=float,
        help="Voxel size (resolution) for the exported mask (nm). (default: %(default)s nm)",
    )
    parser.add_argument(
        "--selection",
        default="not resname W ION",
        type=str,
        help="MDAnalysis selection string for the atom group over which to perform the segmentation. (default: '%(default)s')",
    )
    labels_group = parser.add_mutually_exclusive_group()
    labels_group.add_argument(
        "--labels",
        type=process_labels,
        help="Pre-selected compartment labels, comma-separated.",
    )
    labels_group.add_argument(
        "--autofill",
        action="store_true",
        help="Automatically select the leaf nodes from the containment graph, which commonly represent the 'insides' of the system.",
    )
    parser.add_argument(
        "-b",
        "--inspect-labels-path",
        type=Path,
        nargs="?",
        const="labels.gro",
        help="""Optional output path to a gro file to write labeled voxel positions to.
        This file can be inspected with molecule viewers, which is very helpful 
        in determining which labels match the compartments you want to select.
        (when used, default: %(const)s)""",
    )
    parser.add_argument(
        "--debug-voxels",
        type=Path,
        help="""Write the final voxel mask as a gro file for inspection with molecule viewers. 
        This can be useful when you want to verify the voxel mask that is produced for some selection of labels.""",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="""Do not cache structure files.

        Caching can improve structure load times a lot, because it can load a
        previously stored cache of a MDA Universe. Since loading structure
        files with MDA can be very slow and loading the pickled object is
        relatively quick, this is fantastic for making multiple different
        masks of the same structure.

        This option allows you to *switch off* the caching, if you want that.
        """,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Display verbose output.",
    )
    return parser
