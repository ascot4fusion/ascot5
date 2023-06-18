===========
Data format
===========

Simulation inputs and options, as well as the outputs once the simulation is complete, are stored in a single HDF5 file.

The ASCOT5 HDF5 file, `ascot.h5` from here on out, contains all data required to run a simulation, and so the results and the inputs that produced them are always kept together.
The file is designed to hold multiple inputs (even of same type) and results of multiple simulations to encourage one to use a single file for the whole project or study.

The exact structure of the HDF5 file is irrelevant since it should always be accessed via the Python interface provided by `a5py`.

===========
Data access
===========

Instances of :class:`.Ascot` class are used to open ``ascot.h5``.
The contents of the file are accessed via `Ascot.data`, which is an instance of :class:`.Ascot5IO`, and it provides a "treeview" of the file.
For example, the currently *active* magnetic field input is

.. code-block:: python

   from a5py import Ascot

   a5 = Ascot("ascot.h5")
   a5.data.bfield.active

The magnetic field object here, like the results and input data in general, are stored as objects that inherit from the :class:`.DataGroup` class.
These objects don't contain any data, but provide an interface to access the data on the disk.
This means that initializing the :class:`.Ascot` object is a very light-weight process since no actual data is read during the initialization.

Simulation results are located at the top level of the hierarchy, e.g. the :class:`.DataGroup` corresponding to the *active* run is

.. code-block:: python

   a5.data.active()

Inputs are divided on groups depending on what type of input data they provide. All magnetic field inputs are located under the `bfield` group, all electric field inputs under the `efield`, and so on.
All different *input parent* groups are listed in `INPUGROUPS` and they are instances of :class:`.InputNode` class.
One of the inputs in each group is always flagged *active*, meaning it will be used as an input the next time a simulation is run.
Similarly one of the simulation results is always flagged active (by default it is the latest) and they are instances of :class:`.ResultNode` class.
The active inputs are also used in post-processing and in GUI, so pay attention to what groups are active at the moment.
The group can be set active with

.. code-block:: python

   a5.data.bfield.q0123456789.activate()

The active flag is stored in the HDF5 file so it persists between Python sessions.

Each group has a QID which is a 10-number identifier to separate that group from others and it is automatically created.
For humans there is a description attribute that one can use to document contents of ``ascot.h5``.

The meta data is accessed and modified using the methods found in :class:`.DataGroup` where one can also find methods to remove data groups.
See the tutorial and API below for detailed description of data access.

.. autosummary::
   :nosignatures:

   a5py.ascot5io.Ascot5IO
   a5py.ascot5io.RunGroup
   a5py.ascot5io.InputGroup
   a5py.ascot5io.coreio.treedata.DataGroup

================
Input generation
================

Inputs are written to ``ascot.h5``

.. code-block:: python

   a5 = Ascot("ascot.h5", create=True)

   efield = {"E" : [0,0,0]}
   a5.data.create_input("E_TC", efield)

``E_TC`` is the input type, which in this case is *trivial Cartesian electric field*.
The dictionary provided as the second argument contains data that is passed to an input type specific ``write_hdf5`` function that documents what the data should contain.
The available inputs are listed below.

For input templates and imports, see Importing and Predefined Inputs, but note that some of these use Ascotpy.

.. note::
   ASCOT5 aims to use the user-provided data "as is" as far as possible, so that the user retains control and knowledge what is done with the data.
   This means that ASCOT5 does not extrapolate inputs, so it is user's responsibility to ensure e.g. that the plasma profiles covers the whole region and that the magnetic field extends all the way to the wall.

Magnetic field ``bfield``
=========================

A good quality magnetic field is
If you don't specifically require some other input, use the axisymmetric field since that is fast to interpolate and divergence free.
MHD modes are included via dedicated input and they must not be included in the magnetic field data.

.. autosummary::
   :nosignatures:

   a5py.ascot5io.bfield
   a5py.ascot5io.bfield.B_2DS
   a5py.ascot5io.bfield.B_2DS.write_hdf5
   a5py.ascot5io.bfield.B_3DS
   a5py.ascot5io.bfield.B_3DS.write_hdf5
   a5py.ascot5io.bfield.B_STS
   a5py.ascot5io.bfield.B_STS.write_hdf5
   a5py.ascot5io.bfield.B_GS
   a5py.ascot5io.bfield.B_GS.write_hdf5
   a5py.ascot5io.bfield.B_TC
   a5py.ascot5io.bfield.B_TC.write_hdf5


Electric field ``efield``
=========================

If electric field is not relevant for your simulation, use ``E_TC`` and set it to zero to effectively disable electric field.

.. autosummary::
   :nosignatures:

   a5py.ascot5io.efield
   a5py.ascot5io.efield.E_TC
   a5py.ascot5io.efield.E_TC.write_hdf5
   a5py.ascot5io.efield.E_1DS
   a5py.ascot5io.efield.E_1DS.write_hdf5

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

   a5py.ascot5io.plasma
   a5py.ascot5io.plasma.plasma_1D
   a5py.ascot5io.plasma.plasma_1D.write_hdf5
   a5py.ascot5io.plasma.plasma_1DS
   a5py.ascot5io.plasma.plasma_1DS.write_hdf5

Wall mesh ``wall``
==================

Wall input is required when losses are modelled accurately.

.. autosummary::
   :nosignatures:

   a5py.ascot5io.wall
   a5py.ascot5io.wall.wall_2D
   a5py.ascot5io.wall.wall_2D.write_hdf5
   a5py.ascot5io.wall.wall_3D
   a5py.ascot5io.wall.wall_3D.write_hdf5

Neutral ``neutral``
===================

Neutral particle profiles required when CX reactions are included.

.. autosummary::
   :nosignatures:

   a5py.ascot5io.neutral
   a5py.ascot5io.neutral.N0_3D
   a5py.ascot5io.neutral.N0_3D.write_hdf5

Boozer data ``boozer``
======================

Boozer data is required for simulations with MHD eigenfunctions.
One can create it automatically from ``B_2DS`` with .

.. autosummary::
   :nosignatures:

   a5py.ascot5io.boozer
   a5py.ascot5io.boozer.Boozer
   a5py.ascot5io.boozer.Boozer.write_hdf5

MHD eigenfunctions ``mhd``
==========================

.. autosummary::
   :nosignatures:

   a5py.ascot5io.mhd
   a5py.ascot5io.mhd.MHD_STAT
   a5py.ascot5io.mhd.MHD_STAT.write_hdf5
   a5py.ascot5io.mhd.MHD_NONSTAT
   a5py.ascot5io.mhd.MHD_NONSTAT.write_hdf5

Neutral beam injectors ``nbi``
==============================

.. autosummary::
   :nosignatures:

   a5py.ascot5io.nbi
   a5py.ascot5io.nbi.NBI
   a5py.ascot5io.nbi.NBI.write_hdf5

Markers ``marker``
==================

.. autosummary::
   :nosignatures:

   a5py.ascot5io.marker
   a5py.ascot5io.marker.Prt
   a5py.ascot5io.marker.Prt.write_hdf5
   a5py.ascot5io.marker.GC
   a5py.ascot5io.marker.GC.write_hdf5
   a5py.ascot5io.marker.FL
   a5py.ascot5io.marker.FL.write_hdf5

Options ``opt``
===============

There is only one type of options ``Opt`` and it is treated like any other input.
The options are documented in :ref:`Running Simulations<Simulationoptions>`.

===============
Post-processing
===============

Ascotpy
=======

Many post-processing (and also some pre-processing) tools make use of the Python interface to ASCOT5, referred as ``Ascotpy``.
Instances of ``Ascot`` automatically establish the interface if ``libascot.so`` has been compiled and found in ``LD_LIBRARY_PATH``.
The interface provides direct access to the C routines and it can be used to interpolate inputs exactly as they are interpolated during the simulation.
It can also be used to simulate markers directly from Python.

The interface is used by first initializing the data (note that this can be memory intensive with large inputs).
The `~.Ascot.input_init` method reads the data from the HDF5 file and initializes it using ``ctypes`` so that the data can be passed to C routines.
To free resources, `~.Ascot.input_free` must be called.

.. code-block:: python

   import numpy as np
   from a5py import Ascot

   a5 = Ascot("ascot.h5")
   a5.input_init(bfield=True)
   a5.input_plot(r=np.linspace(4,8,100), z=np.linspace(-2, 2, 100), phi=0, t=0, qnt="rho")
   a5.input_free()

The last line uses the interface internally to evaluate and then plot ``rho``.
An exception is raised if trying to evaluate a quantity that needs input which has not been initialised.

The methods of ``Ascot`` class provide tools to access the input data and they are listed below.

.. autosummary::
   :nosignatures:

   a5py.Ascot
   a5py.Ascot.input_init
   a5py.Ascot.input_free
   a5py.Ascot.input_eval
   a5py.Ascot.input_plotrz

Interactive simulations
=======================

Simulations can be run directly from Python via the interface provided by ``libascot.so`` and ``Ascotpy``.
These simulations are equivalent to running ``ascot5_main`` with the exception that the output is not stored in ``ascot.h5``.
These *interactive* simulations are intented for simulating a handful of markers at most which is why the distribution output is not implemented.
The interface is useful for visualizing marker trajectories or Poincaré plots and evaluating quantities that require markers to be traced only for few orbits.

Before an interactive simulation can be launched, one has to initialize and pack inputs, and set options and markers.

.. code-block:: python

   import numpy as np
   from a5py import Ascot

   a5 = Ascot("ascot.h5")
   a5.simulation_initinputs()
   a5.simulation_initmarkers()
   a5.simulation_initoptions()

The marker and options data are dictionaries with same format as required by their corresponding ``write_hdf5`` routines.
The ``simulation_initinputs`` method requires that no inputs are initialized before calling it.
It differs from the ``input_init`` method in that it initializes all inputs required for the simulation (except options and markers) and *packs* them.
Packing means that input data is allocated in a single array and this is what ASCOT5 does internally when the data is *offloaded* for simulation.
When inputs are packed, the input data cannot be deallocated or altered but otherwise it is possible to use methods like ``input_eval``, that evaluate and plot input data, normally.

Simulation is run and data is accessed as:

.. code-block:: python

   a5.simulation_run()

   # Calling this before the simulation would get the marker inistate
   a5.simulation_getstate()

   # Simulations can be rerun by free'ing the previous results first
   a5.simulation_free(output=True)

   a5.simulation_setoptions()
   a5.simulation_run()
   a5.simulation_getorbit()

   # Free all resources
   a5.simulation_free(inputs=True, markers=True, output=True)

.. autosummary::
   :nosignatures:

   a5py.Ascot
   a5py.Ascot.simulation_initinputs
   a5py.Ascot.simulation_initmarkers
   a5py.Ascot.simulation_initoptions
   a5py.Ascot.simulation_run
   a5py.Ascot.simulation_getstate
   a5py.Ascot.simulation_getorbit
   a5py.Ascot.simulation_free


Simulation output
=================


Ini- and endstate
*****************

Marker orbits
*************

Transport coefficients
**********************

WIP

Distributions
*************

WIP

Losses and wall loads
*********************

End conditions
**************

.. autosummary::
   :nosignatures:

   a5py.ascot5io.state.State.ABORTED
   a5py.ascot5io.state.State.ABORTED
   a5py.ascot5io.state.State.NONE
   a5py.ascot5io.state.State.TLIM
   a5py.ascot5io.state.State.EMIN
   a5py.ascot5io.state.State.THERM
   a5py.ascot5io.state.State.WALL
   a5py.ascot5io.state.State.RHOMIN
   a5py.ascot5io.state.State.RHOMAX
   a5py.ascot5io.state.State.POLMAX
   a5py.ascot5io.state.State.TORMAX
   a5py.ascot5io.state.State.CPUMAX

Exceptions
==========

.. autosummary::
   :nosignatures:

   a5py.exceptions.AscotIOException
   a5py.exceptions.AscotNoDataException
   a5py.exceptions.AscotInitException
