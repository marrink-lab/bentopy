# _hose_&mdash;solvate molecular systems

> **WARNING:** This software is in active development and cannot be considered
> to be reliable or finished. Do not use it expecting usable results.

Solvate systems for molecular dynamics simulations very quickly.

Reads from and writes to [gro files][gro].

## installation

Make sure you have [`rust`][rust] installed on your system. This is very easy. The best method is to use [`rustup`][rustup].

### directly through cargo

```console
$ cargo install --git https://git.sr.ht/~ma3ke/hose
```

### clone this repo, then build

```console
$ git clone https://git.sr.ht/~ma3ke/hose
$ cd hose
$ cargo install --path .
```

## usage

```console
$ hose structure.gro water.gro structure_solvated.gro
```

## Boundary modes

A boundary mode can be set to either cut or grow the final structure box
size.
When the boundary mode is set to `cut`, the original box is kept as it
was provided, and the solvent is cut away in a manner that is aware of
its periodic solvent neigbors. (Default behavior.)
With the boundary mode is set to `grow`, the output structure box is
grown to include the new whole solvent boxes. That means that the output
structure box may expand by at most the size of the template solvent
box, compared to the box set in the input structure.

## Periodicity

The periodic mode determines how the periodic images of the input
structure itself are treated. The default behavior is orthorhombicly
periodic according to the final box size (see boundary mode).
For cases where---for some reason---atoms from the input structure that
fall outside the box should not be considered in checking solvent
placement, the `ignore` mode can be used.
Furher, for cases where solvation should be stopped when an atom from
the input structure exceeds the final boundaries, the boundary mode can
be set to `deny`. This will exit the process with error code 2 if such
an escaped atom is found.

[gro]: https://manual.gromacs.org/archive/5.0.3/online/gro.html
[rust]: https://www.rust-lang.org/
[rustup]: https://www.rust-lang.org/tools/install

---

By Marieke Westendorp, 2024.

<ma3ke.cyber@gmail.com>
