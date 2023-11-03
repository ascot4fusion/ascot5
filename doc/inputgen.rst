.. _Preprocessing:

================
Input generation
================

.. admonition:: Summary

   - Inputs are modular: to simulate markers in a stellarator magnetic field, provide a stellarator magnetic field input.
   - All required inputs must be present even though they would not be actually used: provide dummy atomic data even if you have atomic physics disabled.
   - Inputs are created via :meth:`~.Ascot5IO.create_input`:

     - Dummy input (provide input type but no data): ``create_input("B_2DS")``.
     - Creating input explicitly (provide input type and data): ``create_input("B_2DS", **data)``.
     - Creating input implicitly via templates (provide template name and arguments if needed): ``create_input("bfield analytical iter circular", spline=True)``.

   - Input data can be read with ``read()`` method: ``a5.data.bfield.active.read()``.
   - Information on what inputs there are, what data they need and what templates there are can be found at the end of this page.

ASCOT5 is modular in terms of inputs.
In the code, inputs are interpolated via interfaces so that the implementation can vary.
For example, the magnetic field interface specifies methods to interpolate magnetic field vector at the given coordinates, but the implementation of those methods is different for axisymmetric tokamaks and stellarators.

This modularity is reflected in the format of the HDF5 file, which can look like this for example:

::

    data
    ├── bfield               # Magnetic field inputs
    │   ├── B_2DS_7027705680 # Some axisymmetric tokamak magnetic field data
    │   ├── B_STS_0890178582 # Some stellarator magnetic field data
    │   └── ...              # Some other possible magnetic field data
    │
    ├── efield               # Electric field inputs
    │   └── ...
    │
    └── ...                  # Other inputs (wall, plasma, etc.)

Here ``bfield`` and ``efield`` are parent groups for which there are corresponding magnetic and electric field interfaces in the code.
The actual data is contained in ``B_2DS`` and ``B_STS`` and the code recognizes that these are axisymmetric tokamak and stellarator magnetic field data, respectively, and uses the corresponding implementation in the simulation (depending on which of the inputs was active).

.. note::
   To run a simulation, all required parent groups must have an input of some type present even though the input would not actually be used in the simulation.
   For those cases one can supply dummy data.

Inputs are created with :meth:`~.Ascot5IO.create_input` method which can operate in three ways.
To create *dummy* input, call the method without providing any data.

.. code-block:: python

   # By default the created input is not activated if there are other inputs
   # already present
   a5.data.create_input("E_TC", activate=True)

Here ``E_TC`` is the desired input type, which in this case is the *trivial Cartesian electric field*.

To generate proper input, one must provide all required data which is specified in that input's ``write_hdf5`` function.
In this example, :meth:`.E_TC.write_hdf5` requires that an array with all electric field components is provided.

.. code-block:: python

   a5.data.create_input("E_TC", exyz=np.array([0, 0, 0]))

Existing data can be read with ``read`` method to a dictionary that has the same format as what was passed to ``write_hdf5``:

.. code-block:: python

   data = a5.data.efield.active.read()                   # Read data to a dictionary
   data["exyz"] = exyz=np.array([1, 0, 0])               # Modify
   a5.data.create_input("E_TC", **data, desc="MODIFIED") # Write back to disk

However, it can be laborous to e.g. extract magnetic field data from EQDSK and convert it to the format required by ASCOT5.
For this reason, there are tools called *templates* that convert common data formats to ASCOT5 inputs.
Templates provide a convenient way to import data to ASCOT5 and they also contain some premade inputs, e.g. options and markers for creating Poincaré plots.

Templates are used by providing the name of the template and any required or optional parameters.

.. code-block:: python

   a5.data.create_input("options poincare", maxrho=True, simmode=2, tor=[0, 180])

Next we list all available input types followed by available templates.
When beginning a new study, look through inputs to see what input types you will need and what data is required, and skim through templates to see if you can make use of them.

.. note::
   ASCOT5 aims to use the user-provided data "as is" as far as possible, so that the user retains control and knowledge what is done with the data.
   This means that ASCOT5 does not extrapolate inputs, so it is user's responsibility to ensure e.g. that the plasma profiles covers the whole region and that the magnetic field extends all the way to the wall.

.. currentmodule:: a5py.ascot5io

Magnetic field ``bfield``
=========================

A good quality magnetic field is essential for any orbit-following study.
Always check the field quality by generating Poincaré plots and plotting the divergence.
If you don't specifically require some other input, use the axisymmetric field since that is fast to interpolate and divergence free.
Note that MHD eigenmodes can be included via dedicated input.

.. autosummary::
   :nosignatures:

   bfield
   bfield.B_2DS
   bfield.B_2DS.write_hdf5
   bfield.B_3DS
   bfield.B_3DS.write_hdf5
   bfield.B_STS
   bfield.B_STS.write_hdf5
   bfield.B_GS
   bfield.B_GS.write_hdf5
   bfield.B_TC
   bfield.B_TC.write_hdf5


Electric field ``efield``
=========================

Note that electric field cannot be disabled manually in simulations.
If electric field is not relevant for your simulation, use :class:`.E_TC` and set it to zero to effectively disable electric field.

.. autosummary::
   :nosignatures:

   efield
   efield.E_TC
   efield.E_TC.write_hdf5
   efield.E_1DS
   efield.E_1DS.write_hdf5

..
   These are implemented but not yet merged
   a5py.ascot5io.efield.E_3D
   a5py.ascot5io.efield.E_3D.write_hdf5
   a5py.ascot5io.efield.E_3DS
   a5py.ascot5io.efield.E_3DS.write_hdf5
   a5py.ascot5io.efield.E_3DST
   a5py.ascot5io.efield.E_3DST.write_hdf5

Plasma ``plasma``
=================

Plasma data is required if collisions are included.

.. autosummary::
   :nosignatures:

   plasma
   plasma.plasma_1D
   plasma.plasma_1D.write_hdf5
   plasma.plasma_1DS
   plasma.plasma_1DS.write_hdf5

Wall mesh ``wall``
==================

Wall input is need when modelling losses accurately.
If losses are of no interest it is advised to use some wall to prevent markers from escaping the computational domain.
Alternatively, one can use the ``RHOMAX`` end condition to stop markers at the separatrix.

.. autosummary::
   :nosignatures:

   wall
   wall.wall_2D
   wall.wall_2D.write_hdf5
   wall.wall_3D
   wall.wall_3D.write_hdf5

Neutral ``neutral``
===================

Neutral particle profiles required when atomic (CX, etc.) reactions are included.

.. autosummary::
   :nosignatures:

   neutral
   neutral.N0_1D
   neutral.N0_1D.write_hdf5
   neutral.N0_3D
   neutral.N0_3D.write_hdf5

Boozer data ``boozer``
======================

Boozer data is required for simulations with MHD eigenfunctions.
One can create it automatically from :class:`.B_2DS` or :class:`.B_3DS` with a template :meth:`.InputFactory.boozer_tokamak`.

.. autosummary::
   :nosignatures:

   boozer
   boozer.Boozer
   boozer.Boozer.write_hdf5

MHD eigenfunctions ``mhd``
==========================

MHD input is used to model particle response to MHD (feedback from particles to modes is not implemented).

.. autosummary::
   :nosignatures:

   mhd
   mhd.MHD_STAT
   mhd.MHD_STAT.write_hdf5
   mhd.MHD_NONSTAT
   mhd.MHD_NONSTAT.write_hdf5

Atomic reaction data ``asigma``
===============================

Data for interpolating and computing reaction probabilities and cross sections for reactions where the test particle charge state changes.
This data, typically sourced from ADAS, is required only for simulations where atomic reactions are enabled.

.. autosummary::
   :nosignatures:

   asigma
   asigma.Asigma_loc
   asigma.Asigma_loc.write_hdf5

Neutral beam injectors ``nbi``
==============================

NBI input is used by BBNBI exclusively.
The input consists of a bundle of injectors, which in ``a5py`` are represented by :class:`~nbi.Injector`.

.. autosummary::
   :nosignatures:

   nbi
   nbi.NBI
   nbi.NBI.write_hdf5

Markers ``marker``
==================

Simulation markers.

.. autosummary::
   :nosignatures:

   marker
   marker.Prt
   marker.Prt.write_hdf5
   marker.GC
   marker.GC.write_hdf5
   marker.FL
   marker.FL.write_hdf5

Options ``opt``
===============

There is only one type of options :class:`.Opt` and it is treated like any other input.
The options are documented in :ref:`Running Simulations<Simulationoptions>`.

.. autosummary::
   :nosignatures:

   options
   options.Opt

Templates
=========

.. currentmodule:: a5py.templates

.. rubric:: Data import

.. autosummary::
   :nosignatures:

   ~InputFactory.boozer_tokamak
   ~InputFactory.import_geqdsk
   ~InputFactory.import_adas
   ~InputFactory.mhd_consistent_potentials
   ~InputFactory.ascot4_particles
   ~InputFactory.ascot4_wall2d
   ~InputFactory.ascot4_wall3d
   ~InputFactory.ascot4_tokamak
   ~InputFactory.ascot4_stellarator
   ~InputFactory.ascot4_plasma1d
   ~InputFactory.ascot4_erad
   ~InputFactory.ascot4_neutral1d
   ~InputFactory.ascot4_alfven

.. rubric:: Premade studies

.. autosummary::
   :nosignatures:

   ~InputFactory.options_poincare
   ~InputFactory.marker_poincare
   ~InputFactory.options_tutorial

.. rubric:: Analytic

.. note::
   These are mainly used for testing.

.. autosummary::
   :nosignatures:

   ~InputFactory.bfield_analytical_iter_circular
   ~InputFactory.plasma_flat
   ~InputFactory.neutral_flat
   ~InputFactory.wall_rectangular
