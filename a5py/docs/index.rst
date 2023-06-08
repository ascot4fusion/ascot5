..
  Main page of the ASCOT5 documentation

======
ASCOT5
======


`ASCOT5 <https://arxiv.org/abs/1908.02482/>`_ is a test-particle orbit-following code for solving minority species' distribution functions, transport, and losses in tokamaks and stellarators.
If you can't find here what you are looking for in this documentation, don't be afraid to ask at our `Slack channel <ascot-workspace.slack.com>`_.

.. toctree::
  :maxdepth: 2
  :caption: Running Simulations
  :hidden:

  running
  
.. toctree::
  :maxdepth: 2
  :caption: Processing Data
  :hidden:
  
  processing
  
.. toctree::
  :maxdepth: 2
  :caption: Tutorials
  :hidden:

  teaching
  
.. toctree::
  :maxdepth: 2
  :caption: Publications
  :hidden:

  publishing
  
.. toctree::
  :maxdepth: 2
  :caption: For Developers
  :hidden:

  developing

Getting started
===============

1. Follow the :ref:`installation<Compiling>` instructions and compile the code.
2. Have some quality time by going through the introductory simulation.
3. Check out the examples showcasing different features.
4. Make yourself familiar with the simulation options and input generation to deduce what your work requires.

Features
========

**ASCOT5** is a test-particle orbit-following code for computing particle orbits in 3D geometry. The output includes particle orbits, phase-space distributions, transport coefficients, and wall loads. ASCOT5 is frequently applied to study fast ions, impurities, neutrals, and runaway electrons in tokamaks and stellarators.
Particle orbits are either solved fully, i.e. including the gyro-motion, or in a reduced picture where only the guiding-center trajectory is traced.
The code is extensively parallelized and optimized to support simulations with more than ten million markers.

The code is implemented in C and consists of a main program and a library.
Simulation input and output is stored in HDF5 format.
All pre- and post-processing is done via the Python interface **a5py**. In addition the repository also contain several codes that supplement the orbit-following simulations.

The repository is maintained by Aalto University and VTT.

**ASCOT5:**

- Can trace full gyro-orbits, guiding centers, field-lines, and neutrals.
- Supports 3D tokamak and stellarator magnetic fields.
- Supports 3D wall model and evaluation of particle loads on the wall.
- Simulations assume test-particle approximation (no feedback from markers
  to backgrounds) but otherwise any type of particles can be traced.
  - Mainly used to study fast ions, impurities and runaway electrons.
- Physics include Coulomb collisions and charge-exchange reactions.
- Simulations may include fast ion response to MHD.
- Output consists of various distributions (1D-6D), marker orbits,
  wall loads/FILD signals, and transport coefficients.

**BBNBI5**

- Calculates beam birth profile and shinethrough from NBI geometry.
- Provides NBI source for ASCOT5 slowing-down simulations.

**AFSI**

- Calculates neutron source from fast ion slowing-down distribution
  (as computed by ASCOT5).

**BioSaw**

- Calculates magnetic field based on a coil geometry.
- Provides error fields for ASCOT5 simulations.

**BMC**

- Backward Monte-Carlo simulation tool that is essentially ASCOT5
  but time-reversed.
