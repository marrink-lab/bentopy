import argparse
from pathlib import Path

from render._render import py_render_placements as render_placements

from .grocat import grocat
from .mask import mask
from .check import check
from .pack import pack

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
        type=Path,
        help="""
        Path to the placement list.

        To read from stdin, pass "-".
        """,
    )
    render_parser.add_argument(
        "output",
        type=Path,
        help="Output gro file path.",
    )
    render_parser.add_argument(
        "--root",
        type=Path,
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
    render_parser.add_argument(
        "--resnum-mode",
        type=str,
        choices=["instance", "segment"],
        default="instance",
        help="""Write out a unique resnum for each segment instance, or use one
        grouped resnum for each instance of a segment. (default: %(default)s)""",
    )

    mode_or_topol = render_parser.add_mutually_exclusive_group()
    mode_or_topol.add_argument(
        "-t",
        "--topology",
        dest="topol",
        type=Path,
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

    grocat_parser = subparsers.add_parser(
        "grocat",
        help="Concatenate gro files.",
    )
    grocat_parser.add_argument(
        "files",
        type=grocat.InputFile,
        nargs="+",
        help="""Files to concatenate (gro; <path>[:<resname>]). 

        Optionally, a residue name can be set for all atoms in a file by 
        appending a colon followed by the residue name. 
        Note that this name can be at most 5 characters long. 

        Replacing the residue names can be very useful in distinguishing between 
        parts of very large systems within a concatenated file.""",
    )
    grocat_parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        required=True,
        help="Output path.",
    )
    grocat_parser.add_argument(
        "-t",
        "--title",
        type=str,
        default="bentopy grocat",
        help="Set the final title. (default: %(default)s)",
    )
    grocat_parser.add_argument(
        "-b",
        "--box",
        type=grocat.parse_boxvec,
        help="""Set the final box vectors. 
        Expects a valid gro box line, which is a space-separated list of either 3 or 9 floats. 
        By default, the box vector of the first file is chosen.""",
    )

    check_parser = subparsers.add_parser(
        "check",
        help="Check for collisions in rendered structures.",
    )
    check_parser.add_argument(
        "file",
        type=Path,
        help="File to check (gro)."
    )
    check_parser.add_argument(
        "--cutoff",
        type=float,
        default=0.02,
        help="Collision distance between beads in nm. (default: %(default)s nm)",
    )
    check_parser.add_argument(
        "-e",
        "--exit-early",
        action="store_true",
        help="Exit early at the first collision.",
    )

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
    elif args.subcommand == "grocat":
        grocat.main(args)
    elif args.subcommand == "check":
        check.main(args)
