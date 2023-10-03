..
  Main page of the ASCOT5 documentation

======
ASCOT5
======

https://github.com/ascot4fusion/ascot5

`ASCOT5 <https://arxiv.org/abs/1908.02482/>`_ is a test-particle orbit-following code for solving minority species' distribution functions, transport, and losses in tokamaks and stellarators.
For questions related to the code or physics, please join our `Slack channel <http://ascot-workspace.slack.com/>`_.

.. toctree::
  :maxdepth: 2
  :caption: Running Simulations
  :hidden:

  installing
  simulations

.. toctree::
  :maxdepth: 2
  :caption: Processing Data
  :hidden:

  dataform
  inputgen
  postproc

.. toctree::
  :maxdepth: 2
  :caption: Simulating Physics
  :hidden:

  tutorials
  physics

.. toctree::
  :maxdepth: 2
  :caption: Publications
  :hidden:

  acknowledging
  gallery

.. toctree::
  :caption: Code Development
  :hidden:

  parallel
  testing
  developing

.. toctree::
   :hidden:

   papi/a5py.rst

Getting started
===============

1. Follow the :ref:`installation<Installing>` instructions and :ref:`compile<Compiling>` the code.
2. Have some quality time by going through :ref:`the introductory simulation</tutorials/tutorial.ipynb>`.
3. Familiarize yourself on how to generate :ref:`inputs<Preprocessing>` that you need, execute :ref:`simulations<Simulations>`, and post-process :ref:`the results<Postprocessing>`.
   Here :ref:`the examples<Examples>` and :ref:`the physics<Physics>` documentation as well as :ref:`the Python API</papi/a5py.rst>` are good sources of help.
4. At some point you might also want to :ref:`publish<Citing>` your work or :ref:`contribute<Codingstyle>` to the code.

Features
========

**ASCOT5** is a test-particle orbit-following code for computing particle orbits in 3D geometry.
The output includes particle orbits, phase-space distributions, transport coefficients, and wall loads. ASCOT5 is frequently applied to study fast ions, impurities, neutrals, and runaway electrons in tokamaks and stellarators.
Particle orbits are either solved fully, i.e. including the gyro-motion, or in a reduced picture where only the guiding-center trajectory is traced.
The code is extensively parallelized and optimized to support simulations with more than ten million markers.

The code is implemented in C and consists of a main program and a library.
Simulation input and output is stored in HDF5 format.
All pre- and post-processing is done via the Python interface **a5py**.
In addition the repository also contain several codes that supplement the orbit-following simulations.

The repository is maintained by `Aalto University <https://www.aalto.fi/en/department-of-applied-physics/particle-orbit-simulations-ascot>`_ and `VTT Technical Research Centre of Finland <https://www.vttresearch.com/en/ourservices/future-nuclear#fusion-reactor-research>`_.

**ASCOT5:**

- Can trace full gyro-orbits, guiding centers, field-lines, and neutrals.
- Supports 3D tokamak and stellarator magnetic fields.
- Supports 3D wall model and evaluation of particle loads on the wall.
- Simulations assume test-particle approximation (no feedback from markers to backgrounds) but otherwise any type of particles can be traced.
  - Mainly used to study fast ions, impurities and runaway electrons.
- Physics include Coulomb collisions and charge-exchange reactions.
- Simulations may include fast ion response to MHD.
- Output consists of various distributions (1D-6D), marker orbits, wall loads/FILD signals, and transport coefficients.

**BBNBI5**

- Calculates beam birth profile and shinethrough from NBI geometry.
- Can provides a NBI source for ASCOT5 slowing-down simulations.

**AFSI**

- Calculates fusion product source from thermal plasma and fast ion slowing-down distributions (as computed by ASCOT5).
  For fusion neutronics, this can be combined with `Serpent <https://serpent.vtt.fi/serpent/>`_.

**BioSaw**

- Calculates magnetic field based on a coil geometry.
- Can provide a 3D field for ASCOT5 simulations.

**BMC**

- Backward Monte-Carlo simulation tool that can be thought as a time-reversed ASCOT5.
- Effective tool for estimating FILD signals and wall loads on small but critical components.
- Not yet fully complete.

FAQ and Troubleshooting
=======================

**I'm receiving a fatal error when opening the input file in Python or in GUI**

- There might be compatibility issues if the file was made with an older version of ASCOT5 that you are currently using.
  These can be solved by running the script ``a5update ascot.h5`` which generates a copy ``ascot_<current version>.h5`` where the data has been updated to new format.
  Note that it might not always be possible to enforce the backwards compatibility on results so rerunning the simulation might be necessary.

**How long does it take to run a simulation and how many markers I need**

- See :ref:`this table<ExampleRuntimes>` for examples.
  Typically guiding-center simulations are faster than gyro-orbit because GC can use larger time-step.
  This is somewhat balanced by the fact that a single GC time-step is more expensive to perform.
  Therefore, there is roughly a factor of ten difference between required CPU time but this can vary.
  Using 3D magnetic field instead of 2D is the single most expensive choice in terms of CPU time.
  As for markers, the rule of thumb is that 10k - 100k markers are sufficient to get a good estimate on losses and distributions, but more than 1 million might be needed to estimate wall loads accurately.
  For wall loads the important factors are the size of the wall triangles and how the markers were sampled.
  Always perform some kind of convergence scan or otherwise estimate the Monte Carlo error.
