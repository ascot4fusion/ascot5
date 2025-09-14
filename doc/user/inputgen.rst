.. _Preprocessing:

=================
Generating inputs
=================

The inputs are created using the :class:`~a5py.Ascot.data` class and its various ``create_x`` methods.
`Would you like to know more? <h8.02dasd/>`_

Alternatively, inputs can be created using existing templates, which either:

- Import data (e.g. from EQDSK) and convert it to ASCOT5 format.
- Provide premade inputs for common use cases (e.g. markers and simulation options for Poincaré plots).
- Modify existing inputs (e.g. convert 2D wall to a 3D mesh).

Magnetic field ``bfield``
=========================

.. currentmodule:: a5py.data.bfield

.. toctree::
   :hidden:

   reference/bfieldcartesian
   reference/bfieldanalytical
   reference/bfield2d
   reference/bfield3d
   reference/bfieldstellarator

Magnetic field input is used always in every simulation and it is required to evaluate data via libascot.

A good quality magnetic field is essential for any orbit-following study.
Always check the field quality by generating Poincaré plots and plotting the divergence.
If you don't specifically require some other input, use the axisymmetric field since that is fast to interpolate and divergence free.
Note that MHD eigenmodes can be included via dedicated input.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   BfieldCartesian
   BfieldAnalytical
   Bfield2D
   Bfield3D
   BfieldStellarator


Electric field ``efield``
=========================

.. currentmodule:: a5py.data.efield

.. toctree::
   :hidden:

   reference/efieldcartesian
   reference/efieldradialpotential

Electric field data is used always in every simulation.

Note that electric field cannot be disabled via simulation options.
To do that, generate :class:`.EfieldCartesian` input and set the field to zero.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   EfieldCartesian
   EfieldRadialPotential


Plasma species
==============

``plasma``

.. currentmodule:: a5py.data.plasma

.. toctree::
   :hidden:

   reference/plasma1d
   reference/plasma1ddynamic

Plasma input is used whenever you have ``ENABLE_CCOLL=1``, ``ENABLE_ATOMIC=1``, ``ENDCOND_EMIN=1`` or when you use AFSI.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   Plasma1D
   Plasma1DDynamic


Wall mesh ``wall``
==================

.. currentmodule:: a5py.data.wall

.. toctree::
   :hidden:

   reference/wall2d
   reference/wall3d

Wall input is used when you have ``ENDCOND_WALLHIT=1``.

If losses are of no interest it is advised to use some wall to prevent markers from escaping the computational domain.
Alternatively, one can use the ``RHOMAX`` end condition to stop markers at the separatrix.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   Wall2D
   Wall3D


Neutral species
===============

``neutral``

.. currentmodule:: a5py.data.neutral

.. toctree::
   :hidden:

   reference/neutral1d
   reference/neutral3d

Neutral particle profiles required when atomic (CX, etc.) reactions are included.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   Neutral1D
   Neutral3D


MHD eigenfunctions ``mhd``
==========================

.. currentmodule:: a5py.data.mhd

.. toctree::
   :hidden:

   reference/mhdstationary
   reference/mhddynamic

MHD input is used when ``ENABLE_MHD=1``.

MHD input is used to model particle response to MHD (feedback from particles to modes is not implemented).

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   MhdStationary
   MhdDynamic


Boozer data ``boozer``
======================

.. currentmodule:: a5py.data.boozer

.. toctree::
   :hidden:

   reference/boozer

Boozer input is used when ``ENABLE_MHD=1``.

Boozer data is required for simulations with MHD eigenfunctions.
One can create it automatically from :class:`.B_2DS` or :class:`.B_3DS` with a template :meth:`.InputFactory.boozer_tokamak`.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   BoozerMap


Atomic reaction data ``asigma``
===============================

.. currentmodule:: a5py.data.asigma

.. toctree::
   :hidden:

   reference/atomicdata

Atomic input is used when ``ENABLE_ATOMIC=1``.

Data for interpolating and computing reaction probabilities and cross sections for reactions where the test particle charge state changes.
This data, typically sourced from ADAS, is required only for simulations where atomic reactions are enabled.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   AtomicData


Neutral beam injectors ``nbi``
==============================

.. currentmodule:: a5py.data.nbi

.. toctree::
   :hidden:

   reference/nbi

NBI input is used by BBNBI exclusively.
The input consists of a bundle of injectors, which in ``a5py`` are represented by :class:`~nbi.Injector`.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   NBI


Markers
=======
``marker``

.. currentmodule:: a5py.data.marker

.. toctree::
   :hidden:

   reference/particlemarker
   reference/guidingcentermarker
   reference/fieldlinemarker

Simulation markers are always used except for BBNBI.

.. rubric:: Available inputs

.. autosummary::
   :nosignatures:

   ParticleMarker
   GuidingCenterMarker
   FieldlineMarker


Options
=======
``options``

.. toctree::
   :hidden:

   reference/options

Options are used in every simulation.
There is only one variant of options and it is documented `here <reference/options.rst>`_.

.. tab-set::

   .. tab-item:: Simulation

      .. options-table:: simulation

   .. tab-item:: Physics

      .. options-table:: physics

   .. tab-item:: End conditions

      .. options-table:: endconditions

   .. tab-item:: Distributions

      .. options-table:: distributions

   .. tab-item:: COM-Distribution

      .. options-table:: comdistribution

   .. tab-item:: Orbit

      .. options-table:: orbit

   .. tab-item:: Transport coefficient

      .. options-table:: transport_coefficient


Templates
=========

.. currentmodule:: a5py.templates

Templates are tools that generate ASCOT5 inputs from various sources.
They either use analytical models, import data from external sources (e.g. EQDSK), manipulate existing inputs, or provide premade simulation settings.
We strongly encourage to see if your case is already covered before implementing your own.

.. rubric:: Data import

.. autosummary::
   :nosignatures:

   ~InputFactory.boozer_tokamak
   ~InputFactory.import_geqdsk
   ~InputFactory.import_plasma_profiles
   ~InputFactory.import_adas
   ~InputFactory.import_wall_vtk
   ~InputFactory.import_marsf
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

.. admonition:: Summary

   - Inputs are modular: to simulate markers in a stellarator magnetic field, provide a stellarator magnetic field input.
   - All required inputs must be present even though they would not be actually used: provide dummy atomic data even if you have atomic physics disabled.

     - Information on when input is used can be found at the end of this documentation below each input type's description.

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
