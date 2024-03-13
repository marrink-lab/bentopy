# Learn to pack stuff with _bentopy_

With _bentopy_ you can set up well-stirred systems of biomolecular structures.

![`cool_image.png`]()

## Before we start

First, make sure [you have installed _bentopy_][bentopy_installation].

Note that for this tutorial, it may be convenient to work within a Python virtual environment.

## Contents

- ...

# Introduction

What, why, how, who?

# Preparation of the input files

## Input structures

By design, _bentopy_ is agnostic with respect to force-fields or simulation scale.
The preparation of the structure files is a specific step, in this regard.
All steps following this preparation procedure can be applied to both all-atom and coarse-grained systems.

The goal here, is to prepare the structure files and input topologies such that the packing can be used for a molecular dynamics simulation.

### All-atom

See #9.

### Coarse-grained

Martinize

(A quick EM, eq, md step? Present as optional / irrelevant / do on your own?)


## Input configuration

Make a nice `input.json`. First, step-by-step. Finally, show the finished `input.json`, and we can also just have a correct `input.json` in the `tutorial` directory that you can just download.

# Packing

Demonstrate packing. Explain the `--rearrange` and `--seed` options. Mention the `--verbose` flag, despite its horribly confusing output.

## Inspect the packing

Not implemented, yet. Good place to track the progress of the tools we have for this. (See #12)

# Rendering the placements

Explain what happens in this step.

As an `<aside>`, explain why we do it that way.`</aside>`

## Rendering a small part for inspection

Here be screenshots!

# Energy minimization

Would be good to have all of this set up in a directory within `tutorials`. Just ready to go, ready to take.

## Detecting collisions

(See #2)

# Solvation

Ah, this once again induces me to consider working on a custom way for solvation... (See #11)

## Another EM

# Equilibration 

# A short production run

Here be animations!

[bentopy_installation]: ../README.md#installation
