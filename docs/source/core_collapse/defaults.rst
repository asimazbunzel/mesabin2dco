===========
cc_controls
===========

model_name
~~~~~~~~~~

Name of supernova model to use when a star reaches the core-collapse stage. Possible values are ``rapid``, ``delayed``,
``startrack`` and ``combine``. The first three can be found
`here <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`__ (see Section 4), while a reference for the
latter is `this one <https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.1908K/abstract>`__.

::

    model_name = 'rapid'


max_ns_mass
~~~~~~~~~~~

Maximum value for the mass of a neutron-star in solar masses. Needed for the prescriptions found in Fryer et al. paper.

::

    max_ns_mass = 2.5d0


report_core_collapse
~~~~~~~~~~~~~~~~~~~~

Flag to report on the terminal information of the core-collapse calculations

::

    report_core_collapse = .false.


cc_data_directory
~~~~~~~~~~~~~~~~~

Name of the direcotry where core-collapse information will be saved.

::

    cc_data_directory = 'cc_data'


filename_for_star_data
~~~~~~~~~~~~~~~~~~~~~~

Name of the file containing information of the collapsing star at the core-collapse stage.

::

    filename_for_star_data = 'star_at_core_collapse'


filename_for_binary_data
~~~~~~~~~~~~~~~~~~~~~~~~

Name of the file with information of the binary system at core-collapse.

::

    filename_for_binary_data = 'binary_at_core_collapse'
