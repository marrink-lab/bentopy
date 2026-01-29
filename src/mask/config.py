from typing import Literal
import argparse
from pathlib import Path

EPILOG = """
Compartments are determined using the MDVContainment package created by Bart Bruininks.
More information is available at <https://github.com/BartBruininks/mdvcontainment>.
"""
DESCRIPTION = "Set up masks based on structures and compartment segmentations."

DEFAULT_CONTAINMENT_RESOLUTION = 0.5
DEFAULT_MASK_RESOLUTION = 0.5


def parse_label(label: str) -> int | None:
    match label:
        case "autofill": return None
        case n: return int(n)


def process_labels_to_mask(labels_to_mask: str) -> tuple[list[int | None], str]:
    """
    Convert a string of the form '<label>+:<mask-path>' to
    `(labels, mask_path)`.

    For example, '-1,-2:mask.npz' becomes `([-1, -2], "mask.npz")`.

    A possible label is 'autofill', which will be interpreted during runtime
    based on the containment graph.

    For example, 'autofill,-2:mask.npz' becomes `([None, -2], "mask.npz")`.
    """
    labels_list, mask_path = labels_to_mask.split(":")
    labels = [parse_label(label.strip()) for label in labels_list.split(",")]
    return labels, mask_path


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
        "-l",
        "--labels",
        nargs=1, # Make sure that we can also safely provide negative number arguments.
        action="append",
        type=process_labels_to_mask,
        help="""Write compartments to some mask based on a list of labels.
        Provide compartment labels, comma-separated, followed by a colon and
        the mask path for that set of labels.

        For example, providing '-1,-2:mask.npz' will create a mask called
        'mask.npz' that represents the labeled compartments -1 and -2.
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
    parser.add_argument(
        "-b",
        "--visualize-labels",
        type=Path,
        help="""Optional output path to a gro file to write labeled voxel
        positions to. This file can be inspected with molecule viewers, which
        is very helpful in determining which labels match the compartments you
        want to select.""",
    )
    parser.add_argument(
        "--exclude-outside",
        action="store_true",
        help="""When writing labeled voxel positions (--visualize-labels),
        exclude empty outside nodes.""",
    )
    parser.add_argument(
        "--visualize-masks",
        action="store_true",
        help="""Write the final voxel masks as a gro file for inspection with
        molecule viewers. This can be useful when you want to verify the voxel
        mask that is produced for some selection of labels.

        The name will be derived from the output mask name. So 'mask.npz'
        becomes 'mask.npz.gro'.
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
    morph = parser.add_mutually_exclusive_group()
    morph.add_argument(
        "--morph",
        default="de",
        type=str,
        help="""Morphological operations to apply to the initial boolean voxel
        representation based on the provided structure. Provide a string of 'd'
        for dilation steps and 'e' for erosion steps, in the order you wish to
        apply them.

        The default 'de' is equivalent to a morphological closing operation,
        while 'ed' is a morphological opening operation. See
        <https://en.wikipedia.org/wiki/Closing_(morphology)> for more
        information.

        If the resolution exceeds the condensed phase distance (i.e., about
        double the LJ sigma), morphing is not required. A voxel resolution
        below sigma is not recommended.
        (default: '%(default)s')
        """,
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

def parse_args():
    parser = setup_parser()
    args = parser.parse_args()
    # Because we use both nargs=1 and action="append", we need to unpack our
    # values. Sucks but it works.
    if args.labels is not None:
        args.labels = [l for [l] in args.labels]

    # Make sure that --exclude-outside is accepted iff --visualize-labels is
    # given.
    if args.exclude_outside and not args.visualize_labels:
        parser.error("--exclude-outside requires --visualize-labels")

    return args
