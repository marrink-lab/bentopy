# _bentopy_&mdash;packs stuff in boxes

## State

This project is in an early state. Many parts are in flux, and no guarantees
about correctness or stability can be made.

## Installation

### Prerequisites

As of now, _bentopy_ requires a [Rust][rust] compiler to be installed. To check
whether this is the case, you can run

```console
cargo --version
```

If it is not present, [you can install it][rust-installation] by any means you
prefer. Generally, installation through _rustup_ is most convenient.

### Install _bentopy_ through _pip_ directly

If you don't care about peeking into the sources and just want access to the
program, this is the quickest option.

```console
python3 -m venv venv && source venv/bin/activate # Not required, but often convenient.
pip3 install git+https://github.com/marrink-lab/bentopy
```

### From source

```console
git clone https://github.com/marrink-lab/bentopy
cd bentopy
python3 -m venv venv && source venv/bin/activate
pip3 install .
```

## Usage

_bentopy_ currently features two subcommands, _pack_ and _render_.

You can learn about the available options through the help information.

```console
bentopy --help
```

### _pack_

_pack_ provides the core functionality of _bentopy_. Given an input
configuration file, a packing of the input structures is created. 

This packing is stored as a **placement list**, which is a `json` file that
describes which structures at what rotations are placed where. In order to
create a structure file from this placement list that can be read by molecular
visualazation and simulation programs, the _render_ subcommand can be used.

---

Let's try to pack a system that is full of lysozyme structures. First, we want
to get a structure to pack. Let's download the structure for [`3LYZ`][3lyz] and
put it into a `structures` directory to stay organized.

```console
wget https://files.rcsb.org/download/3lyz.pdb
mkdir structures
mv 3lyz.pdb structures
```

Now we can set up our input configuration, which we will call `3lyz_input.json`:

```json
{
	"space": {
		"size": [100, 100, 100],
		"resolution": 0.5,
		"compartments": [
			{
				"id": "main",
				"shape": "spherical"
			}
		]
	},
	"output": {
		"title": "3lyz",
		"dir": "output",
		"topol_includes": [
			"forcefields/forcefield.itp",
			"structures/3lyz.itp"
		]
	},
	"segments": [
		{
			"name": "3lyz",
			"number": 6500,
			"path": "structures/3lyz.pdb",
			"compartments": ["main"]
		}
	]
}
```

We set the **space** up to a **size** of 100&times;100&times;100 nm, with a
**resulution** of 0.5 nm. The mask&mdash;the volume that defines where
structures can be placed&mdash;is set to be derived from a **spherical** analytical function.

In case you want to use a custom **mask**, you could specify one in the following manner.

```json
		"compartments": [
			{
				"id": "main",
				"voxels": {
					"path": "inputs/sphere.npz"
				}
			}
		]
```

Here, **voxels** and the associated **path** point to a precomputed voxel mask.
This mask can be any data that can be loaded by [`np.load()`][np-load] to be
interpreted as a three-dimensional boolean mask. The provided mask must have
the same shape as specified in the **space** section's **dimensions**.
By specifying custom masks, a great range of systems can be packed!

> [!CAUTION]
> Currently, multiple masks (beyond the single arbitrarily named "main" from
> this example) are either broken or very unstable.
> See #7.

In **output**, we set a **title** and **dir**ectory to write the product files
to. With the optional field **topol_includes**, we can specify what 
[`itp` files][itp] files are to be included if the placement list produced from
this config is written to a topology file. 

> [!NOTE]
> For this example, we do not actually fill this field. It is shown here for
> the purpose of explanation, not as a useful example.
> See #9.

Finally, in the **segments** section, we define a list of structures to place.
In our case that is only one: which we give the **name** "3lyz", and we set the
**number** of segments to place to 6500. The **path** points _pack_ to where
the structure file for this segment can be found.

> [!IMPORTANT]
> The **name** record must be selected carefully. If you want to write out a
> valid topology file using _render_, the value of **name** must correspond to
> the names in the `itp` files.

Now, we are ready to pack the system. Most simply, we can do this as follows.

```console
bentopy pack 3lyz_input.json
```

In order to make the procedure deterministic, the `--seed` parameter can be set.
This means that the same command will produce the same output between runs.

```console
bentopy pack --seed 1312 3lyz_input.json
```

After the command finishes, we will find that `output/3lyz_placements.json` has
been created. This is a single-line `json` file, which can be hard to inspect.
If you are curious, you can use a tool such as [`jq`][jq] to look at what was
written.

```console
jq . output/3lyz_placement.json
```

<details>
<summary>
The output may look like this (some lines have been cut and adjusted for legibility).
</summary>

```
{
	"title": "3lyz",
	"size": [ 100, 100, 100 ],
	"topol_includes": [ ... ],
	"placements": [
		{
			"name": "3lyz",
			"path": "structures/3lyz.pdb",
			"batches": [
				[
					[
						[ 1.0, 0.0, 0.0 ], 
                        [ 0.0, 1.0, 0.0 ], 
                        [ 0.0, 0.0, 1.0 ]
					],
					[
						[  8, 46, 68 ],
						[ 26, 62, 88 ],
                        ... many many more of such lines ...
                    ]
                ],
				[
					[
						[   0.3658391780537972, -0.3882572475566672, -0.8458238619952991  ],
						[  -0.8851693094147572, -0.4258733932991502, -0.18736901171236636 ],
						[ -0.28746650147647396,  0.8172442490465064, -0.49947457185455224 ]
					],
					[
						[ 31, 41, 56 ],
						[ 61, 53,  4 ],
                        ... many many more of such lines ...
                    ]
                ]
                ... and on and on and on ...
            ]
        }
    ]
}
```
</details>

### _render_

_render_ reads in the placement list and writes out a [`gro` file][gro] 
(and optionally, a [`top` topology file][top]). This is a separate operation,
since the packed systems can become very large. Storing the placement list as
an intermediate product decouples the hard task of packing from the simple work
of writing it into a structure file.

---

We want to render out the placement list we just created into a structure file
called `3lyz_sphere.gro`. Additionally, we would like to produce topology file
(`topol.top`) that Gromacs uses to understand how the structure file is built
up.

```console
bentopy render output/3lyz_placements.json 3lyz_sphere.gro -t topol.top
```

You can now inspect the `3lyz_sphere.gro` structure in a molecular visualization
program of your preference.

But beware! We just created big structure, and some programs may have a hard
time keeping up. Luckily, _bentopy render_ has some additional tricks up its
sleeve to ease this load.

In case you want to inspect only a small part of a very large placement list,
the `--limits` option allows you to select a cuboid within the volume defined
by the placement list from which the placed structures will be rendered. The
volume that is cut out is defined by a sequence of six comma-separated values
in the order `minx,maxx,miny,maxy,minz,maxz`. If a value is a number, it is
interpreted as a dimension in nm. If it is not a number (the phrase 'none' is
conventional) no limits are set on that dimension. 

For example, to only render a 10&times;10&times;10 nm cube extending
from the point (40, 40, 40) to (50, 50, 50), we can pass the following limits.

```console
bentopy render output/3lyz_placements.json 3lyz_small_cube.gro --limits 40,50,40,50,40,50
```

Perhaps we would like to see a pancake instead! To do this, we can define the
limits only for the _z_-direction.

```console
bentopy render output/3lyz_placements.json 3lyz_pancake.gro --limits none,none,none,none,45,55
```

Using `--limits`, we can cut out a part of the packed structure, but perhaps
you want to inspect the total structure without loading as many atoms.

For this, you can try the `--mode` option, which gives you the ability to only
render out certain atoms (`backbone`, `alpha` carbon) or beads (representing
each `residue`, or even only one per structure `instance`). By default, the
mode is `full`, and we have just seen its output. Let's try `alpha`, now.

> [!WARNING] 
> Some of these options (`backbone`) may not be functional right now.

```console
bentopy render output/3lyz_placements.json 3lyz_alpha.gro --mode alpha
```

Now, we can compare the sizes of the files.

```console
wc -l 3lyz_sphere.gro 3lyz_alpha.gro
```

Reducing the number of atoms that are rendered out can improve the time it
takes to inspect a packing, if necessary.

> [!NOTE]
> Using modes other than `full` (the default) is obviously not relevant beyond
> inspection and analysis of the packed structure. To reflect this, the option
> to write a topology file and setting a mode are mutually exclusive.

In case you want to render out a structure based on a placement list that you
or a colleague have created in a different environment, it can be useful to
direct _render_ to read the input structures from a different directory. To do
this, you can set a root path for the structures. This path will be prepended
to any relative structure path that is defined in the placement list.

[rust]: https://rust-lang.org/
[rust-installation]: https://www.rust-lang.org/learn/get-started
[gro]: https://manual.gromacs.org/current/reference-manual/file-formats.html#gro
[top]: https://manual.gromacs.org/current/reference-manual/file-formats.html#top
[top]: https://manual.gromacs.org/current/reference-manual/file-formats.html#itp
[3lyz]: https://www.rcsb.org/structure/3LYZ
[jq]: https://github.com/jqlang/jq
[np-load]: https://numpy.org/doc/stable/reference/generated/numpy.load.html
