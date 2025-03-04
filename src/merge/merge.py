import argparse
from sys import stderr

# This should be sufficient space to write the final natoms into.
NATOMS_PLACEHOLDER = " " * 32
RESNAME_MAX_LEN = 5


class InputFile:
    divider = ":"

    def __init__(self, s: str):
        if self.divider in s:
            path, resname = s.split(self.divider)
            if len(resname) > RESNAME_MAX_LEN:
                raise ValueError(
                    f"A resname cannot be longer than 5 characters, found '{resname}' with {len(resname)} characters"
                )
        else:
            path = s
            resname = None
        self.file = argparse.FileType("r")(path)
        self.resname = resname


def eprint(*args, **kwargs):
    """
    Print to standard error.
    """
    print(*args, **kwargs, file=stderr)


def guarantee_newline(line: str) -> str:
    """
    If the provided line does not have a trailing newline, yet, add it.
    Otherwise leave it alone.

    Note that this is an idempotent operation (i.e., ∀x f(x) = f(f(x))).
    """
    if line.endswith("\n"):
        return line
    else:
        return line + "\n"


def parse_boxvec(s: str) -> tuple[float]:
    """
    Parse a gro box vector line.
    """
    boxvec = [float(v) for v in s.split()]
    nv = len(boxvec)
    if nv == 3 or nv == 9:
        return boxvec
    raise ValueError(
        f"Could not parse '{s}'. The box vector must have either 3 or 9 values, found {nv}."
    )


def merge(args):
    output = args.output
    if not output.seekable():
        eprint("ERROR: Output into non-seekable files is currently not supported.")
        return 1

    # Write the title.
    if args.title is not None:
        output.write(guarantee_newline(args.title))

    boxvecs = []
    natoms_total = 0
    natoms_location = None  # The location in output where the natoms are written.
    # We know that there will be at least one file provided through args.
    for input_file in args.files:
        # TODO: Unpack this with a tuple in the for loop? How can we implement that for a custom class?
        file = input_file.file
        resname = input_file.resname

        # Report on the file and if we replace the resnames.
        if resname is None:
            info = "Leaving resnames in place."
        else:
            assert (
                len(resname) <= RESNAME_MAX_LEN
            )  # Already verified during argument parsing.
            info = f"Replacing resnames with '{resname}'."
        eprint(f"Reading from {file.name}...", info)

        next(file)  # Skip the title.

        # Put in a natoms placeholder if this is the first file that is written.
        if natoms_location is None:
            natoms_location = output.tell()
            output.write(NATOMS_PLACEHOLDER + "\n")

        # Read the number of atoms that are in this file.
        natoms = int(next(file).strip())

        # Write that number of lines into the output file.
        if resname is not None:
            # Modify the resname for each atom.
            padded_resname = f"{resname:<5}"
            for _ in range(natoms):
                line = next(file)
                output.write(line[:5])
                output.write(padded_resname)
                output.write(line[10:])
        else:
            for _ in range(natoms):
                output.write(next(file))

        # Now that we have written the lines, we'll add them to our total.
        natoms_total += natoms

        # Parse the box vector and add it to our collection.
        boxvec = parse_boxvec(next(file))
        boxvecs.append(boxvec)

    # TODO: Take into account args.box.
    # Write out the final box vector.
    final_boxvec = boxvecs[0]
    output.write(" ".join([str(v) for v in final_boxvec]) + "\n")

    # Finally, we go and seek back to the start, where we left our natoms placeholder.
    # Let's replace that with the correct total number of atom records we wrote to output.
    output.seek(natoms_location)
    output.write(str(natoms_total))


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate gro files.",
        prog="bentopy-merge",
    )

    parser.add_argument(
        "files",
        type=InputFile,
        nargs="+",
        help="""Files to concatenate (gro; <path>[:<resname>]). 

        Optionally, a residue name can be set for all atoms in a file by 
        appending a colon followed by the residue name. 
        Note that this name can be at most 5 characters long. 

        Replacing the residue names can be very useful in distinguishing between 
        parts of very large systems within a concatenated file.""",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        required=True,
        help="Output path.",
    )
    parser.add_argument(
        "-t",
        "--title",
        type=str,
        default="bentopy merge",
        help="Set the final title. (default: %(default)s)",
    )
    parser.add_argument(
        "-b",
        "--box",
        type=parse_boxvec,
        help="""Set the final box vectors. 
        Expects a valid gro box line, which is a space-separated list of either 3 or 9 floats. 
        By default, the box vector of the first file is chosen.""",
    )

    args = parser.parse_args()
    merge(args)


if __name__ == "__main__":
    main()
