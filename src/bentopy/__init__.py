import argparse

from extensions._extensions import py_voxelize as voxelize

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

    args = parser.parse_args()

    if args.version:
        import importlib.metadata

        version = importlib.metadata.version("bentopy")
        print(version)
        return 0

    if args.subcommand == "pack":
        state = pack.configure(args)
        pack.main(state)
