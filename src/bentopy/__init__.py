import argparse

from extensions._extensions import py_voxelize as voxelize

from .grocat import grocat
from .pack import pack

__all__ = ["voxelize"]


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

    args = parser.parse_args()

    if args.version:
        import importlib.metadata

        version = importlib.metadata.version("bentopy")
        print(version)
        return 0

    if args.subcommand == "pack":
        state = pack.configure(args)
        pack.main(state)
    elif args.subcommand == "grocat":
        grocat.main(args)
