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

- Now is time to compile the whole source code. Do that by running `./mk`.

- Finally, to run the simulation just type `./rn`.

## Options and controls

### bin2dco

### Common-envelope phase

#### alpha_two_stars

Efficiency for the removal of envelope of a donor star in a binary containing two non-degenerate stars.

#### alpha_xrb

Efficiency for the removal of envelope of a donor star in a binary containing a non-degenerate star and a compact object.

#### edd_scaling_factor

Eddington scaling factor used to trigger the common-envelope phase.

#### tol_two_stars

Tolerance to consider common-envelope detachment when having two non-degenerate stars. A binary is considered to detach
after the common-envelope phase if `(r-rl)/rl` is below this number. Here `r` is the radius of the donor star and `rl`
the corresponding Roche lobe.

#### tol_xrb

Tolerance to consider common-envelope detachment when having a star and a compact object. Same definition as for the
`tol_two_stars` applies.

#### years_to_max_mdot_rlof

Increment the value of `mdot_rlof` during this many years until reach `max_mdot_rlof`.

#### max_mdot_rlof

Maximum value for the mass-transfer through the Roche lobe overflow.

#### add_accretion_on_ce

Flag to consider accretion during common-envelope. **NOT READY TO USE**.

#### save_profile_pre_ce
#### save_profile_after_ce

Flags to control the storage of MESA profiles just before and after (if possible) the common-envelope phase.

#### save_model_pre_ce
#### save_model_after_ce

Flags to control the storage of MESA models just before and after (if possible) the common-envelope phase.

#### ce_data_directory

Name of the folder where common-envelope information will be saved.

#### filename_donor_profile_pre_ce
#### filename_donor_profile_after_ce

Names of files for the MESA profiles of the donor star just before and after (if possible) the common-envelope phase.

#### filename_donor_model_pre_ce
#### filename_donor_model_after_ce

Names of files for the MESA models of the donor star just before and after (if possible) the common-envelope phase.

#### filename_accretor_profile_pre_ce
#### filename_accretor_profile_after_ce

Names of files for the MESA profiles of the accretor star just before and after (if possible) the common-envelope phase.
Only used when binary consists of two non-degenerate stars.

#### filename_accretor_model_pre_ce
#### filename_accretor_model_after_ce

Names of files for the MESA models of the accretor star just before and after (if possible) the common-envelope phase.
Only used when binary consists of two non-degenerate stars.

#### max_relative_gap

If `(r-rl)/rl` is bigger than this number, then end simulation and consider it a common-envelope merger.

#### max_number_retries_during_ce

If the number of retries exceeds this limit, then end simulation and consider it a common-envelope merger.

### Core-collapse

#### model_name

Name of supernova model to use when a star reaches the core-collapse stage. Possible values are `rapid`, `delayed`,
`startrack` and `combine`. The first three can be found [here](https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract)
(see Section 4), while a reference for the latter is [this](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.1908K/abstract).

#### cc_data_directory

Name of the folder where core-collapse information will be saved.

#### filename_for_star_data

Name of the file containing information of the collapsing star at the core-collapse stage.

#### filename_for_binary_data

Name of the file with information of the binary system at core-collapse.

#### max_ns_mass

Maximum value for the mass of a neutron-star in solar masses. Needed for the prescriptions found in Fryer et al. paper.

#### add_asymmetric_kick

Add a single random kick from inside the core-collapse module. **NOT READY TO USE**.

#### continue_binary_evolution

Flag to control the continuation of a simulation after a binary reaches the core-collapse stage. **NOT READY TO USE**.
