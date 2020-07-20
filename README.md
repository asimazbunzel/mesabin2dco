# mesa-bin2dco
---

## Description

`mesa_bin2dco` evolves two stars, or a star plus a point-mass companion, all the way through the main-sequence and up
to the formation of a double-compact object system.

## Walkthrough of a run

## Requirements

It requires having the open-access, evolutionary code MESA (see their [official website](http://mesa.sourceforge.net/))
installed. In particular, this code is thought to be applied with revision 10398, together with the software-development
kit version 20180822 (get it [here](http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk-old#linux-download)).

## Getting started

- If this is the first time running the code, modules stored in the `src` folder need to be compiled. In order to do so,
  run `./mk_mods` which is found in the root of the tree structure.

- Next, change MESA parameters in the different inlist files as convenient, and also have a look and change controls in
  the `bin2dco_controls` file.

- Now is time to compile the whole source code. Do that by running `./mk`

- Finally, to run the simulation just type `./rn`.

## Options and controls

### bin2dco

### Common-envelope phase

### Core-collapse
