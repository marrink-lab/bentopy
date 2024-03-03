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
$ hose --center structure.gro water.gro structure_solvated.gro
```

[gro]: https://manual.gromacs.org/archive/5.0.3/online/gro.html
[rust]: https://www.rust-lang.org/
[rustup]: https://www.rust-lang.org/tools/install

---

By Marieke Westendorp, 2024.

<ma3ke.cyber@gmail.com>
