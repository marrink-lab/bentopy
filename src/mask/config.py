import argparse
from pathlib import Path

EPILOG = """
Compartments are determined using the MDVContainment package created by Bart Bruininks.
More information is available at <https://github.com/BartBruininks/mdvcontainment>.
"""
DESCRIPTION = "Set up masks based on structures and compartment segmentations."

DEFAULT_CONTAINMENT_RESOLUTION = 0.5
DEFAULT_MASK_RESOLUTION = 0.5


def process_labels(labels: str) -> list[int]:
    return [int(label) for label in labels.split(",")]


def setup_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
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
    morph = parser.add_mutually_exclusive_group()
    morph.add_argument(
        "--morph",
        type=str,
        help="""Morphological operations to apply to the initial boolean voxel
        representation based on the provided structure. Provide a string of 'd'
        for dilation steps and 'e' for erosion steps, in the order you wish to
        apply them.

        For example, the string 'de' is equivalent to the --closing flag, while
        'ed' is a morphological opening operation. See
        <https://en.wikipedia.org/wiki/Closing_(morphology)> for more
        information.
        """,
    )
    morph.add_argument(
        "--closing",
        action="store_true",
        help="""Use binary closing to fill small holes in compartments. For
        analysis of CG Martini structures with a voxel resolution 0.5 nm this
        is highly recommended.

        If the resolution exceeds the condensed phase distance (i.e., about
        double the LJ sigma), closing is not required. A voxel resolution below
        sigma is not recommended.

        See also, --morph, which provides a more flexible interface to the same
        notion.
        """,
    )
    parser.add_argument(
        "--min-size",
        type=float,
        help="""
        A volume in nm³. Starting from the leaf nodes, recursively merge with
        ancestors until this minimum size has been met.

        This is very helpful for filtering out small, non-relevant compartments
        by merging them to their parent compartments.
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
