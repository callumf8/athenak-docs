AthenaK documentation
=====================

.. image:: /_static/athenak_hero.png
   :alt: AthenaK simulation
   :align: center
   :width: 100%

AthenaK is a block-based AMR framework with fluid, particle, and numerical relativity
solvers implemented with `Kokkos <https://kokkos.org/>`_ for performance portability.

This site covers how to obtain, build, and run the code, and documents the
AMR framework and Physics solvers.

.. toctree::
   :caption: Getting Started
   :maxdepth: 1

   requirements
   download
   build

.. toctree::
   :caption: Quickstart
   :maxdepth: 1

   quickstart

.. toctree::
   :caption: Running
   :maxdepth: 1

   running_the_code
   input_file
   outputs
   analysis
   notes_specific_machines

.. toctree::
   :caption: Grid and Initial Conditions
   :maxdepth: 1

   uniform_mesh
   smr
   amr
   cyclic_zoom
   problem_generators

.. toctree::
   :caption: Physics Modules
   :maxdepth: 1

   hydro
   mhd
   ion_neutral
   diffusion
   source_terms
   radiation
   nr_radiation
   numerical_relativity
   dyngrmhd
   dyngrmhd_eos
   particles

.. toctree::
   :caption: Other Features
   :maxdepth: 1

   shearing_box
   units
   autotest

.. toctree::
   :caption: Help
   :maxdepth: 1

   faq
