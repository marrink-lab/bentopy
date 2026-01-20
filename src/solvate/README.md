# _bentopy solvate_&mdash;solvate molecular systems

> Previously, _bentopy solvate_ was known as _hose_. It used to be developed at
> [git.sr.ht/~ma3ke/hose][hose].

Solvate systems for molecular dynamics simulations very quickly, including ion
placement and charge neutralization.

This program is designed to solvate extremely large systems while retaining a
low memory footprint and without wasting time.

Reads from and writes to [gro files][gro].

## Installation

This is a subcommand of _bentopy_, and is installed along with it. Please refer
to this project's [installation instructions][install].

## Usage

```console
$ bentopy-solvate -i structure.gro -o structure_solvated.gro
```

Learn about the available options through `$ bentopy-solvate --help`.

```
Solvate

Usage: bentopy-solvate [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>
          Structure input path

  -o, --output <OUTPUT>
          Output path

      --cutoff <CUTOFF>
          Lowest allowable distance between solvent and structure beads (nm).

          This is the minimum allowable center-to-center distance when checking
          for a collision between a solvent bead and a structure bead.

          For the solvent-solvent cutoff distance, see `--solvent-cutoff`.

          For Martini solvation, the default cutoff is 0.43 nm. For atomistic
          solvation, the default cutoff is 0.28.

      --solvent-cutoff <SOLVENT_CUTOFF>
          Lowest allowable distance between solvent beads (nm).

          This is the minimum allowable center-to-center distance between
          solvent beads when cutting the sides to fit in the specified output
          structure box.

          For Martini solvation, the default cutoff is 0.21 nm. For atomistic
          solvation, the default cutoff is 0.21.

      --ignore <IGNORE>
          List of resnames to ignore when checking against structure-solvent
          collisions

      --water-type <WATER_TYPE>
          The type of water written to the output file

          [default: martini]
          [possible values: martini, tip3p]

  -c, --center
          Center the structure in the new box.

          Note that this option only has an effect if the boundary mode is set
          to `grow`.

  -b, --boundary-mode <BOUNDARY_MODE>
          Set how the boundaries of the box will be treated

          Possible values:
          - cut:  Cut at the structure box size and remove solvent residues
                  that overlap with the periodic neighbors
          - grow: If necessary, grow the box of the output structure to fit a
                  whole number of template boxes

          [default: cut]

  -p, --periodic-mode <PERIODIC_MODE>
          Set how periodic images of the input structure are considered.

          Note that only orthorhombic (cuboid) periodicity is considered,
          currently.

          Possible values:
          - periodic: Include the periodic images of the structure when
                      checking whether a solvent spot is occupied
          - ignore:   Ignore the periodic images of the structure for the
                      solvent placement check
          - deny:     Treat any structure atoms outside of the output box as an
                      error and exit

          [default: periodic]

  -s, --substitute <SUBSTITUTES>
          Define substitutes for solvent residues, such as ions

      --charge <CHARGE>
          Set the charge to neutralize with additional ions.

          In essence, this is a shorthand for explicitly providing ions to
          compensate the charge as substitutes. This can be helpful for
          automation purposes.

          By default, NA (positive) and CL (negative) ions are used, but a
          different positive-negative substitution pair can be specified by
          prepending a colon followed by the positive and negative substitute
          name separated by one comma.

          <charge>:<positive ion name>,<negative ion name>

          So, by default, some `<charge>` is interpreted as

          <charge>:NA,CL

      --no-combine-substitutes
          Combine substitutes with identical names into one block

      --sort-substitutes <SORT_SUBSTITUTES>
          Set whether and in what way substitutes are sorted

          [default: size]
          [possible values: size, rev-size, alphabetical, rev-alphabetical, no]

      --seed <SEED>
          Random number generator seed for solvent substitution, such as ion
          placement

      --write-velocities
          If the solvent template contains velocity information, write these
          velocities to the output file

  -t, --append-topol <APPEND_TOPOL>
          Append solvation topology lines to a path

      --no-write-parallel


      --buffer-size <BUFFER_SIZE>
          The suggested number of atoms to format at once.

          Setting this to a larger value will allow more atoms to be formatted
          at once, at the cost of higher memory consumption. Setting this to a
          smaller value will lower the memory footprint at a possible cost of
          efficiency.

          Note that depending on the number of residues in the template box,
          the provided value may be overshot by at most that number.

          [default: 10000000]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Features

### All-atom and Martini solvation

The water model used for solvation can be selected with the `--water-type`
flag. The `martini` and `tip3p` (atomistic) options are available. By default,
the `martini` water type is used.

### Substitutions&mdash;place ions with ease!

Distributing ions in solvated systems can become challenging in very large
systems. By substituting solvent residues with ions or other residues in this
solvation process, a great speed-up is possible over other methods.

```console
$ bentopy-solvate -i structure.gro -o solvated.gro -s NA:0.15M -s CL:0.15M
```

To specify a substitution, provide a name and a quantifier, separated by a
colon: `<name>:<quantity>`.

A quantity can be expressed in any of the following ways:

- `<float>M`: _molarity_, determine the number of substitutions based on the
  volume of the output structure dimensions.
- `<float>`: _ratio_, the fraction of valid solvent beads to replace.
- `<integer>`: _number_, a set number of substitutions to make.

#### Neutralizing structure charge

The `--charge` option can be used as a convenient short-hand for adding
neutralizing ions. This may be especially useful in automated pipelines.

```console
$ bentopy-solvate -i structure.gro -o solvated.gro -s NA:0.15M -s CL:0.15M --charge -12
```

`NA` and `CL` substitutes are made by default, but custom neutralizing
substitutes can be set after the charge:
`--charge <charge>:<positive name>,<negative name>`.

```console
$ bentopy-solvate -i structure.gro -o solvated.gro -s NA:0.15M -s CL:0.15M --charge -12:K,CL
```

#### Combining, sorting

By default, the substitute records with identical names are combined and sorted
by decreasing size. This combining behavior can be turned off with the
`--no-combine-substitutes` flag, and the sorting behavior can be set with the
`--sort-substitutes` option.

### Boundary modes

A boundary mode can be set to either cut or grow the final structure box
size.
When the boundary mode is set to `cut`, the original box is kept as it
was provided, and the solvent is cut away in a manner that is aware of
its periodic solvent neigbors. (Default behavior.)
With the boundary mode is set to `grow`, the output structure box is
grown to include the new whole solvent boxes. That means that the output
structure box may expand by at most the size of the template solvent
box, compared to the box set in the input structure.

### Periodicity

The periodic mode determines how the periodic images of the input
structure itself are treated. The default behavior is orthorhombicly
periodic according to the final box size (see boundary mode).
For cases where&mdash;-for some reason&mdash;-atoms from the input structure that
fall outside the box should not be considered in checking solvent
placement, the `ignore` mode can be used.
Furher, for cases where solvation should be stopped when an atom from
the input structure exceeds the final boundaries, the boundary mode can
be set to `deny`. This will exit the process with error code 2 if such
an escaped atom is found.

[hose]: https://git.sr.ht/~ma3ke/hose
[install]: /README.md#installation
[gro]: https://manual.gromacs.org/archive/5.0.3/online/gro.html
[rust]: https://www.rust-lang.org/
[rustup]: https://www.rust-lang.org/tools/install

---

By Marieke Westendorp, 2024.

<ma3ke.cyber@gmail.com>
