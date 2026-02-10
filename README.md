# _bentopy_&mdash;packs stuff in boxes

![Bentopy](figures/logo_header.png)

Bentopy packs molecules to assemble models for molecular dynamics simulations.
It uses a voxel-based approach to place molecules into arbitrary geometries,
handling systems from simple boxes to cellular-scale models. Additional tools
provide efficient solvation and ion placement. Bentopy is force-field agnostic
and integrates with existing molecular dynamics workflows. Example systems
built with bentopy include entire cell models in the Martini force field and
all-atom aerosols.

- Build models at vast scales, quickly.
- Powerful and flexible integration of spatial information.
- Works with any forcefield, from Martini to all-atom.

## Information

A number of resources explaining how to build systems using bentopy are
available.

- The [_bentopy_ wiki][wiki] gives in-depth descriptions of different commands
  and concepts.
	- This includes [Examples][wiki-examples] of how different _bentopy_
	  tools can be used to construct simple and more sophisticated systems.
- The in-depth tutorial: [Bentopy: from simple packing to building cellular
  models][tutorial]. This is an updated version of the tutorial presented at the
  [Martini Workshop 2025][workshop].

[wiki]: https://github.com/marrink-lab/bentopy/wiki
[wiki-examples]: https://github.com/marrink-lab/bentopy/wiki/Examples
[tutorial]: https://cgmartini.nl/docs/tutorials/Martini3/Bentopy
[workshop]: https://cgmartini.nl/docs/tutorials/Martini3/workshop

## Citation

If you use bentopy to set up your molecular dynamics system in a publication,
please cite our work.

> **Westendorp MSS**, **Stevens JA**, Brown CM, Dommer AC, Wassenaar TA,
> Bruininks BMH, and Marrink SJ. Compartment-guided assembly of large-scale
> molecular models with _Bentopy_. _Protein Science_. 2026. e70480.
> <https://doi.org/10.1002/pro.70480>

```bibtex
@article{westendorp_stevens_bentopy_2026,
	title = {Compartment-Guided Assembly of Large-Scale Molecular Models with Bentopy},
	doi = {10.1002/pro.70480},
	pages = {e70480},
	journaltitle = {Protein Science},
	author = {Westendorp, Marieke S.S. and
                  Stevens, Jan A. and
                  Brown, Chelsea M. and
                  Dommer, Abigail C. and
                  Wassenaar, Tsjerk A. and
                  Bruininks, Bart M.H. and
                  Marrink, Siewert J.},
	year = 2026,
}
```

## Installation

_Bentopy_ can be installed through `pip`.

```console
pip install bentopy==1.0.0b6
```

> [!NOTE]
> As you can see, we recommend you install a beta version. This will soon
> become a _bentopy_ 1.0 release candidate, and eventually the proper 1.0
> release.
>
> This is the version described by the current documentation.

For most Linux platforms, pre-built binaries are available and will be
installed automatically through `pip`. For other platforms, such as macOS, a
Rust compiler is required. See the [detailed installation
instructions](#detailed-installation-instructions) below.

## The _bentopy_ tools

_Bentopy_ currently features five subcommands:

- [_init_](#init): Initialize and validate bentopy input files.
- [_mask_](#mask): Create masks based on voxel containments.
- [_pack_](#pack): Pack a space and produce a placement list.
- [_render_](#render): Write a structure file and topology based on a placement list.
- [_solvate_](#solvate): Solvate large models, including very fast ion substitution.
- [_merge_](#merge): Merge structure files.

Each subcommand has detailed `--help` information.

### Workflow

A typical _bentopy_ workflow may look like this.

![bentopy-merge -> bentopy-mask -> bentopy-pack -> bentopy-render -> bentopy-merge -> bentopy-solvate](figures/bentopy_workflow.png)

The _init_ and _mask_ tools help prepare the configuration and masks for
packing. Based on the input configuration and masks, _pack_ will create a
placement list. Using _render_, you can create a structure file and system
topology from the placement list. The _merge_ and _solvate_ tools help you
prepare the model for simulation.

What follows is a brief explanation and example invocation of these
subcommands. A more detailed walkthrough can be found in the
[Examples][wiki-examples] on the wiki.

### Pre-processing

#### _init_

The _init_ subcommand serves to make setting up a new _bentopy_ project easy.
It can be used to create an [example configuration file][example] with
placeholder values or to validate input files.

[example]: https://github.com/marrink-lab/bentopy/blob/main/src/init/example.bent

```console
bentopy-init example -o input.bent
```

Read more about the [_init_ command](https://github.com/marrink-lab/bentopy/wiki/bentopy-init).

#### _mask_

While simple shapes can be defined directly in a `bent` file, **space masks**
enable you to capture the complex geometries of curved membranes and large
complexes. This allows you to create models based on existing structures or
empirical evidence, making sophisticated integrative modeling workflows
possible.

The _mask_ subcommand offers a powerful tool for creating these masks. It is
built on top of a versatile library for segmenting point clouds and molecular
structures, called [mdvcontainment][mdvc].

With _mask_ you can take a structure or point cloud and determine the different
compartments within it.

```console
bentopy-mask vesicle.gro -l autofill:masks/inside.npz
```

_Determine the compartments contained by the structure in `vesicle.gro` and
automatically select the innermost compartment (`autofill`). From that
selected compartment, write a mask to `masks/inside.npz`._

The masks created with `bentopy-mask` can be imported as a compartment in a
`bent` file.

```ini
[ compartment ]
cytoplasm from "masks/inside.npz"
```

Note that any boolean _numpy_ array [stored as a compressed file
(`npz`)][numpy-npz] of the correct dimensions can function as a valid mask.
This makes it possible to create custom scripts and techniques for preparing
masks as well.

### Packing the structure

#### _pack_

The _pack_ subcommand provides the core functionality of _bentopy_. Given an
**input configuration file** (`bent`), the input structures will be packed and
their positions and orientations are written to a **placement list**.

```console
bentopy-pack input.bent placements.json
```

The _placement list_ can be converted to a structure and associated topology
using [_render_](#render).

#### The `bent` input configuration file

_On the wiki, a [detailed reference] for the `bent` configuration input file is
available._

[detailed reference]: https://github.com/marrink-lab/bentopy/wiki/Reference-for-bent

This is a minimal but complete `bent` file.

```ini
[ general ]
title "Lysozymes in a vesicle"

[ space ]
dimensions 50, 50, 50
resolution 0.5

[ compartments ]
inside from "masks/inside.npz"

[ segments ]
3lyz 300 from "structures/3lyz.pdb" in inside
```

Provided you have the structure file and the mask, this configuration can be
used to place lysozyme structures within a vesicle.

For a more elaborate description of a similar system using an analytical
sphere, take a look at [Example 1: Simple sphere][wiki-example-1] on the wiki.

[wiki-example-1]: https://github.com/marrink-lab/bentopy/wiki/Example-1:-Simple-sphere

### Post-processing

#### _render_

The result of the packing process is stored as a **placement list**, which is a
`json` file that describes _which structures_ at _what rotations_ are _placed
where_. In order to create a structure file (and topology file) from this
placement list, the _render_ subcommand can be used.

```console
bentopy-render placements.json packed.gro -t topol.top
```

_Render `placements.json` created by _pack_ to a `gro` file at `packed.gro`
and write a topology file to `topol.top`._

This is a separate operation from _packing_, since the packed systems can
become very large. Storing the placement list as an intermediate result
decouples the hard task of packing from the simple work of writing it into a
structure file.

#### _merge_

As the name suggests, _merge_ is a tool for concatenating `gro` files. Though
this is a relatively simple operation, _merge_ provides a convenient way of
telling apart different sections of large models by optionally specifying a new
residue name for a whole file in the argument list by appending
`:<residue name>` to a file path.

```console
bentopy-merge vesicle.gro:VES packed.gro:LYZ -o merged.gro
```

_Concatenate `vesicle.gro` and `packed.gro` into `merged.gro`, setting the
residue names of the lipid vesicle atoms to `VES` and those of the packed
lysozymes to `LYZ` in the concatenated structure._

#### _solvate_

With _solvate_, large boxes can be solvated quickly and conveniently. _Solvate_
places waters and makes ion substitutions in the same step. It was created to
enable cellular-scale solvation and is designed to run very fast while having a
low memory footprint. Both atomistic and coarse-grained Martini water placement
is supported.

```console
bentopy-solvate -i merged.gro -o solvated.gro -t topol.top \
	-s NA:0.15M -s CL:0.15M --charge neutral
```

_Solvate the structure in `packed.gro` with Martini waters and output the
result to `solvated.gro`. Substitute water residues for ions at 0.15M NaCl.
Make additional ion substitutions to neutralize the charge of `packed.gro`
based on the provided topology. The number of added waters and ions will be
appended to `topol.top`._

A thorough description of the command is [provided in the `bentopy-solvate`
README](src/solvate/README.md).

## Detailed installation instructions

If pre-built binaries are not available for your platform, you need access to
`cargo`. First, you can check if Rust's build system `cargo` is installed.

```console
cargo --version
```

Make sure that this is at least `1.88`, the minimum supported rust version for
this project. If this is not the case, update using `rustup update`.

If `cargo` is not present, [you can install it][rust-installation] by any means
you prefer. Installation through [_rustup_][rust-rustup] is very easy!

Once `cargo` is installed, installing _bentopy_ using `pip` should work.

### Install from source

Installing bentopy from source gives you access to the very latest changes.
Though the main branch of this project is generally stable, installation from
releases is recommended.

To install _bentopy_ from source, you need access to cargo, as described
earlier. You can use `pip` directly to install the project right from the
repository.

```console
pip install git+https://github.com/marrink-lab/bentopy
```

Alternatively, you can clone the repository somewhere, and build from there.

```console
git clone https://github.com/marrink-lab/bentopy
cd bentopy
pip install .
```

## License

> Copyright 2024 Marieke S.S. Westendorp, Jan A. Stevens
> 
> Licensed under the Apache License, Version 2.0 (the "License");
> you may not use this file except in compliance with the License.
> You may obtain a copy of the License at
> 
> <http://www.apache.org/licenses/LICENSE-2.0>
> 
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS,
> WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
> See the License for the specific language governing permissions and
> limitations under the License.

[rust]: https://rust-lang.org/
[rust-installation]: https://www.rust-lang.org/learn/get-started
[rust-rustup]: https://rustup.rs/
[numpy-npz]: https://numpy.org/doc/stable/reference/generated/numpy.savez.html
[mdvc]: https://github.com/BartBruininks/mdvcontainment
[gromacs-gro]: https://manual.gromacs.org/current/reference-manual/file-formats.html#gro
[gromacs-top]: https://manual.gromacs.org/current/reference-manual/file-formats.html#top
[gromacs-itp]: https://manual.gromacs.org/current/reference-manual/file-formats.html#itp
[3lyz]: https://www.rcsb.org/structure/3LYZ
[1ubq]: https://www.rcsb.org/structure/1UBQ
[jq]: https://github.com/jqlang/jq
[np-load]: https://numpy.org/doc/stable/reference/generated/numpy.load.html
