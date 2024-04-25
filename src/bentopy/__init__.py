import argparse
import pathlib

from render._render import py_render_placements as render_placements
from .pack import pack
from .mask import mask

__all__ = ["render_placements"]


def main():
    parser = argparse.ArgumentParser(
        description="Packs stuff in boxes.",
        prog="bentopy",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="store_true",
        help="Print the version and exit",
    )

    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    pack_parser = subparsers.add_parser(
        "pack",
        help="""
        Pack a space with structures and produce a placement list of this packing.
        """,
    )
    pack.setup_parser(pack_parser)

    render_parser = subparsers.add_parser(
        "render",
        help="""
        Render structures from a placement list into a gro file.

        Structures specified in the placement list are retrieved from their pdb
        files and placed into a gro file according to their rotations and positions.
        """,
    )
    render_parser.add_argument(
        "input",
        type=pathlib.Path,
        help="""
        Path to the placement list.

        To read from stdin, pass "-".
        """,
    )
    render_parser.add_argument(
        "output",
        type=pathlib.Path,
        help="Output gro file path.",
    )
    render_parser.add_argument(
        "--root",
        type=pathlib.Path,
        help="""
        Root path for the structure paths.
        
        When set, this path will be prepended to any relative path pointing to a structure in the
        placement list. Absolute paths are respected.
        """,
    )
    render_parser.add_argument(
        "--limits",
        type=str,
        help="""
        Only render structures that have a position within a smaller cuboid.

        Arguments can be provided as a comma-separated array of 6 values. Each value can be a
        number indicating a bound or a non-numerical value indicating an unset bound.

        For example, `none,none,none,none,90,110` produces a 20 nm slice.
        """,
    )
    mode_or_topol = render_parser.add_mutually_exclusive_group()
    mode_or_topol.add_argument(
        "-t",
        "--topology",
        dest="topol",
        type=pathlib.Path,
        help="Write a topology (.top) file.",
    )
    mode_or_topol.add_argument(
        "--mode",
        type=str,
        choices=["full", "backbone", "alpha", "residue", "instance"],
        help="Granularity of the produced output.",
    )

    mask_parser = subparsers.add_parser(
        "mask",
        help=mask.DESCRIPTION,
    )
    mask.setup_parser(mask_parser)

    args = parser.parse_args()

    if args.version:
        import importlib.metadata
        import sys

        version = importlib.metadata.version("bentopy")
        print(version)
        sys.exit(0)

    if args.subcommand == "pack":
        state = pack.configure(args)
        pack.main(state)
    elif args.subcommand == "render":
        render_placements(
            args.input, args.output, args.topol, args.root, args.limits, args.mode
        )
    elif args.subcommand == "mask":
        mask.main(args)
