===========
ce_controls
===========

alpha_two_stars
~~~~~~~~~~~~~~~

Efficiency for the removal of envelope of a donor star in a binary containing two non-degenerate stars.

::

    alpha_two_stars = 1d0


alpha_xrb
~~~~~~~~~

Efficiency for the removal of envelope of a donor star in a binary containing a non-degenerate star and a compact object.

::

    alpha_xrb = 2d0


edd_scaling_factor
~~~~~~~~~~~~~~~~~~

Eddington scaling factor used to trigger the common-envelope phase.

::

    edd_scaling_factor = 1d0


tol_two_stars
~~~~~~~~~~~~~

Tolerance to consider common-envelope detachment when having two non-degenerate stars. A binary is considered to detach
after the common-envelope phase if `(r-rl)/rl` is below this number. Here `r` is the radius of the donor star and `rl`
the corresponding Roche lobe.

::

    tol_two_stars = -1d-2


tol_xrb
~~~~~~~

Tolerance to consider common-envelope detachment when having a star and a compact object. Same definition as for the
`tol_two_stars` applies.

::

    tol_xrb = -1d-1


years_in_detachment
~~~~~~~~~~~~~~~~~~~

Binary is detached for these many years before reaching common-envelope end

::

    years_in_detachment = 1d0


years_to_max_mdot_rlof
~~~~~~~~~~~~~~~~~~~~~~

Increment the value of `mdot_rlof` during this many years until reach `max_mdot_rlof`.

::

    years_to_max_mdot_rlof = 7.5d0


max_mdot_rlof
~~~~~~~~~~~~~

Maximum value for the mass-transfer through the Roche lobe overflow.

::

    max_mdot_rlof = 0.1d0


add_accretion_on_ce
~~~~~~~~~~~~~~~~~~~

Flag to consider accretion during common-envelope. *NOT READY TO USE*.

::

    add_accretion_on_ce = .false.


save_profile_pre_ce
~~~~~~~~~~~~~~~~~~~

save_profile_after_ce
~~~~~~~~~~~~~~~~~~~~~

Flags to control the storage of MESA profiles just before and after (if possible) the common-envelope phase.

::

    save_profile_pre_ce = .true.
    save_profile_after_ce = .true.

save_model_pre_ce
~~~~~~~~~~~~~~~~~

save_model_after_ce
~~~~~~~~~~~~~~~~~~~

Flags to control the storage of MESA models just before and after (if possible) the common-envelope phase.

::

    save_model_pre_ce = .true.
    save_model_after_ce = .true.


ce_data_directory
~~~~~~~~~~~~~~~~~

Name of the folder where common-envelope information will be saved.

::

    ce_data_directory = 'ce_data'

filename_donor_profile_pre_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename_donor_profile_after_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Names of files for the MESA profiles of the donor star just before and after (if possible) the common-envelope phase.

::

    filename_donor_profile_pre_ce = 'profile_donor_pre_ce'
    filename_donor_profile_after_ce = 'profile_donor_after_ce'


filename_donor_model_pre_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename_donor_model_after_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Names of files for the MESA models of the donor star just before and after (if possible) the common-envelope phase.

::

    filename_donor_model_pre_ce = 'donor_pre_ce'
    filename_donor_model_after_ce = 'donor_after_ce'

filename_accretor_profile_pre_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename_accretor_profile_after_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Names of files for the MESA profiles of the accretor star just before and after (if possible) the common-envelope phase.
Only used when binary consists of two non-degenerate stars.

::

    filename_accretor_profile_pre_ce = 'profile_accretor_pre_ce'
    filename_accretor_profile_after_ce = 'profile_accretor_after_ce'

filename_accretor_model_pre_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename_accretor_model_after_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Names of files for the MESA models of the accretor star just before and after (if possible) the common-envelope phase.
Only used when binary consists of two non-degenerate stars.

::

    filename_accretor_model_pre_ce = 'accretor_pre_ce'
    filename_accretor_model_after_ce = 'accretor_after_ce'

max_relative_gap
~~~~~~~~~~~~~~~~

If `(r-rl)/rl` is bigger than this number, then end simulation and consider it a common-envelope merger.

::

    max_relative_gap = 1d2

max_number_retries_during_ce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the number of retries exceeds this limit, then end simulation and consider it a common-envelope merger.

::

    max_number_retries_during_ce = 200
