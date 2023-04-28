=======
Modules
=======

Documentation of modules developed for the *mesabin2dco* wrapper.

The wrapper itself is also explain as if it was a module, even though formally it's not a MESA
module. On the other hand, both ``ce`` and ``core_collapse`` are proper MESA modules as they follow
the general structure described by the MESA developers in
`Module documentation <https://docs.mesastar.org/en/release-r15140/modules.html>`__

Bin2dco (``bin2dco``)
---------------------

The ``bin2dco`` wrapper handles the launch of one/many binary stellar-evolution(s).

.. toctree::
   :maxdepth: 1

   bin2dco/overview
   bin2dco/defaults


Core-collapse (``core_collapse``)
---------------------------------

The ``core_collapse`` module handles the computation of the compact object in a collapsing star

.. toctree::
   :maxdepth: 1

   core_collapse/overview
   core_collapse/defaults


Common-envelope (``ce``)
------------------------

The ``ce`` module handles the evolution of a binary during a common-envelope phase

.. toctree::
   :maxdepth: 1

   ce/overview
   ce/defaults
