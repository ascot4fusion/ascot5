.. _Preprocessing:

================
Input generation
================

The inputs in ASCOT5 are all stored in ``ascot.h5`` and its structure reflects how the code operates.
Each input belongs to one of the *parent groups* that in the code itself act like interfaces.
For example, the magnetic field parent group ``bfield`` has a corresponding interface in the code that provides routines to evaluate magnetic field quantities at the queried coordinates etc.
The exact implementation is left for the modules that correspond to the *data groups* located under the parent group in ``ascot.h5``.
This allows ASCOT5 to have modular inputs so that we can have e.g. one implementation for an axisymmetric tokamak field, :class:`.B_2DS`, and one implementation for stellarators, :class:`.B_STS`.
To run a simulation, all required input parent groups must have a data group of some type present even though the input would not actually be used in the simulation.
For those cases one can supply dummy data.

Inputs are created with :meth:`~.Ascot5IO.create_input` method which can operate in three ways.
To create *dummy* input, call the method without providing any data.

.. code-block:: python

   a5.data.create_input("E_TC", desc="Dummy data", activate=True)

Here ``E_TC`` is the desired input type, which in this case is the *trivial Cartesian electric field*.
Description is user-specified, where the first word acts as a *tag* when traversing the treeview (here the tag becomes DUMMY as covered in the last section).
The last argument sets the input active.

To generate proper input, one must provide all required data.
**What data is required can be learned by inspecting the input specific** ``write_hdf5`` **function.**
In this example, :meth:`.E_TC.write_hdf5` requires that an array with all electric field components is provided.

.. code-block:: python

   a5.data.create_input("E_TC", exyz=np.array([0, 0, 0]))

However, it can be laborous to e.g. extract magnetic field data from EQDSK and convert it to the format required by the axisymmetric tokamak field :meth:`.B_2DS.write_hdf5`.
For this reason, we aim to implement conversion of common data formats to ASCOT5 inputs as *templates* that process the data into a suitable format and calls proper ``write_hdf5`` internally.
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

MHD eigenfunctions ``asigma``
=============================

Atomic reaction data (e.g. from ADAS) for simulations where atomic reactions are enabled.

.. autosummary::
   :nosignatures:

   asigma
   asigma.Asigma_loc
   asigma.Asigma_loc.write_hdf5

Neutral beam injectors ``nbi``
==============================

NBI input is used by BBNBI5 exclusively.
The input consists of a bundle of injectors, which in ``a5py`` are represented by :class:`Injector`.

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

.. autosummary::
   :nosignatures:

   ~InputFactory.bfield_analytical_iter_circular
   ~InputFactory.boozer_tokamak
   ~InputFactory.options_poincare
   ~InputFactory.marker_poincare
   ~InputFactory.plasma_flat
   ~InputFactory.wall_rectangular
   ~InputFactory.options_tutorial
