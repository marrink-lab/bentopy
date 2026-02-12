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
          Define substitutes for solvent residues, such as ions.

          A substitute is defined according to the scheme <name>:<quantity>.
          For example, 150 mM NaCl can be described as
            `-s NA:0.15Ms -s CL:0.15Ms`.

          A shorthand for this is `-s NA,CL:0.15Ms`. Note that this shorthand
          respects stoichiometry for ions such as MgCl₂: `-s MG,CL@2:0.026Ms`. The
          `@2` here is used to indicate the stoichiometric ratio of Mg:Cl::1:2 for
          the dissociated magnesium chloride.

          Quantities can be specified as follows.

          - Molar concentration with respect to solvent quantity:
            floating point number followed by an 'Ms' suffix.
            (Example: `-s NA:0.15Ms` replaces 150 mM of solvent residues with
            NA, determined based on the remaining solvent quantity.)

          - Molar concentration with respect to box volume: floating point
            number followed by an 'M' suffix.
            (Example: `-s NA:0.15M` replaces 150 mM of solvent residues with NA,
            determined based on the box volume.) This behaviour is equivalent to
            that of many other solvation tools, such as `gmx genion`.

          - Count: an unsigned integer.
            (Example: `-s NA:100` replaces 100 solvent residues with NA.)

          - Ratio: a floating point number.
            (Example: `-s NA:0.1` replaces 10% of solvent residues are replaced
            with NA.)

      --charge <CHARGE>
          Neutralize the system charge with additional ions.

          By passing `--charge neutral`, the total system charge will be
          determined automically. This requires a topology.

          The charge to be neutralized can also be set explicitly by providing
          an integer.

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

          Here, `<charge>` can be the string 'neutral' or an integer.

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
$ bentopy-solvate -i structure.gro -o solvated.gro -s NA:0.15Ms -s CL:0.15Ms
```

To specify a substitution, provide a name and a quantifier, separated by a
colon: `<name>:<quantity>`.

#### Substitute names

The name can be a single residue name, or a set of residue names separated by
commas, like `-s NA,CL:0.15Ms`. In order to describe stoichiometric ratios of
salts, you can use the `@<number>` notation. For example, a 0.026 M MgCl₂ can be
described as `-s MG,CL@2:0.026Ms`. The `@2` here is used to indicate the
stoichiometric ratio of Mg:Cl::1:2 for the dissociated magnesium chloride.

#### Quantity

A quantity can be expressed in any of the following ways:

- `<float>Ms`: _solvent molarity_, determine the number of substitutions based
  remaining solvent quantity.

  $$n_\text{X} = \frac{c_\text{X} \cdot n_\text{solvent}}{c_\text{X} + c_\text{water}}$$
  where $c_\text{water} = 55.34 \ \text{M}$ at 1 bar, 300 K.
- `<float>M`: _box volume molarity_, determine the number of substitutions based
  on the volume of the output structure dimensions. Note that this is
  classically how many solvation/ion placement tools interpret molarity, such as
  [`gmx genion`][genion]. Especially for systems with a large non-solvent volume
  (say, a membrane with an embedded protein), computing the ion counts based on
  the box volume may lead to a much higher effective ion concentration than
  expected.

  $$n_\text{X} = c_\text{X} \cdot V_\text{box} \cdot N_\text{A}$$
- `<float>`: _ratio_, the fraction of valid solvent beads to replace.
- `<integer>`: _number_, a set number of substitutions to make.

[genion]: https://gitlab.com/gromacs/gromacs/-/blob/release-2026/src/gromacs/gmxpreprocess/genion.cpp#L553

#### Neutralizing structure charge

By passing `--charge neutral`, you can automatically neutralize the total
system charge when a topology is provided.

```console
$ bentopy-solvate -i structure.gro -o solvated.gro -t topol.top \
        -s NA:0.15M -s CL:0.15M --charge neutral
```

The `--charge` option can also be used as a convenient short-hand for adding
neutralizing ions for an explicitly provided integer charge. This may be
especially useful in automated pipelines.

```console
$ bentopy-solvate -i structure.gro -o solvated.gro \
        -s NA:0.15M -s CL:0.15M --charge -12
```

`NA` and `CL` substitutes are made by default, but custom neutralizing
substitutes can be set after the charge:
`--charge <charge>:<positive name>,<negative name>`.

```console
$ bentopy-solvate -i structure.gro -o solvated.gro \
        -s NA:0.15M -s CL:0.15M --charge neutral:K,CL
$ bentopy-solvate -i structure.gro -o solvated.gro \
        -s NA:0.15M -s CL:0.15M --charge -12:K,CL
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
