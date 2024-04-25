from sys import exit, stderr

# This should be sufficient space to write the final natoms into.
NATOMS_PLACEHOLDER = " " * 32


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


def main(args):
    output = args.output
    if not output.seekable():
        eprint("ERROR: Output into non-seekable files is currently not supported.")
        exit(1)

    title = args.title
    # Write the title if we already know what it is.
    if title is not None:
        output.write(guarantee_newline(title))

    boxvecs = []
    natoms_total = 0
    natoms_location = None  # The location in output where the natoms are written.
    # We know that there will be at least one file provided through args.
    for file in args.files:
        eprint(f"Reading from {file.name}...")

        # Set the title if it has not been set, yet.
        if title is None:
            title = next(file)
            output.write(title)
        else:
            next(file)  # Skip this title since we already set it.

        # Put in a natoms placeholder if this is the first file that is written.
        if natoms_location is None:
            natoms_location = output.tell()
            output.write(NATOMS_PLACEHOLDER + "\n")

        # Read the number of atoms that are in this file.
        natoms = int(next(file).strip())
        # Write that number of lines into the output file.
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
