==========
Quickstart
==========

``mesabin2dco`` is a wrapper of the `MESA <https://docs.mesastar.org>`__ that handles the evolution
of

- Two non-degenerate stars (*star_plus_star*) this kind of evolution is defined by the
  ``MESAbinary`` module with the binary_job option ``evolve_both_stars = .true.``.

- Star plus a compact companion (*star_plus_pm*): the code follows the evolution of only one star
  in a close binary with a compact object (``evolve_both_stars = .false.``).

In addition, this wrapper contains two modules that handles different evolutionary stages during
binary evolution: a common envelope phase and a core-collapse stage, ``ce`` and ``core_collapse``
modules respectively.

Installation
------------

In order for this wrapper to work, we need MESA installed, version *r15140*. See their
`Installing MESA <https://docs.mesastar.org/en/release-r22.11.1/installation.html>`__ section.
On top of that, the `MESASDK <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`__
version compatible (tested) with this wrapper is *20.12.1*.

Once this is installed, running ``./mk`` in the root directory of ``mesabin2dco`` will create the
executable ``binary`` that will launch the stellar evolution simulation.

Getting started
---------------

An example is given, with all the ``inlist*`` files, similar to what MESA does. The only difference
is that controls specific for the ``mesabin2dco`` are located in the ``inlist`` file, under the
fortran namelist ``bin2dco_controls``. To see a detail explanation of all these controls, see
:ref:`Controls of mesabin2dco <bin2dco_controls_defaults>`

In addition, there are two other controls to set that are specific for the ``ce`` and
``core_collapse`` modules. Their explanations are located in:
:ref:`Controls of ce <ce_controls_defaults>` and
:ref:`Controls of core_collapse <core_collapse_defaults>`

To run the code there are three different cases which can be used, depending on the type of
evolution:

1. `star_plus_star` and `star_plus_pm`: simply run ``./rn``

2. `star_plus_pm` with a single asymmetric kick from the kth row a file: ``./rn k``

.. note::

   Option 1. will also work in the case of running several asymmetric kicks from a file (in loop)
