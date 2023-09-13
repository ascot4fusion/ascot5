===========
Data format
===========

Simulation inputs and options, as well as the outputs once the simulation is complete, are stored in a single HDF5 file.

The ASCOT5 HDF5 file, ``ascot.h5`` from here on out, contains all data required to run a simulation, and so the results and the inputs that produced them are always kept together.
The file is designed to hold multiple inputs (even of same type) and results of multiple simulations to encourage one to use a single file for the whole project or study.

The exact structure of the HDF5 file is not that relevent since it should always be accessed via the Python interface provided by :mod:`a5py`.

===========
Data access
===========

Instances of :class:`.Ascot` class are used to open ``ascot.h5`` or create an empty file.

.. code-block:: python

   from a5py import Ascot
   a5 = Ascot("ascot.h5", create=True) # Creates a new file

The contents of the file are accessed via ``Ascot.data``, which is an instance of :class:`.Ascot5IO`, and it provides a "treeview" of the file.
For example, the currently *active* magnetic field input is

.. code-block:: python

   a5 = Ascot("ascot.h5")
   a5.data.bfield.active # Currently active magnetic field input

The magnetic field object here, like the results and input data in general, are stored as objects that inherit from the :class:`.DataGroup` class.
These objects don't actually contain any data.
Instead, they provide an interface to access the data on the disk.
This makes :class:`.Ascot` objects very light-weight since no actual data is read during the initialization.

Simulation results are located at the top level of the hierarchy, e.g. the :class:`.DataGroup` corresponding to the *active* run is

.. code-block:: python

   a5.data.active # Currently active run (most recent by default)

Inputs are divided to groups depending on what type of input data they provide.
All magnetic field inputs are located in the ``bfield`` group, all electric field inputs in the ``efield`` group, and so on.
All different *input parent* groups are instances of :class:`.InputNode` class.
**One of the inputs in each group is always flagged *active*, meaning it will be used as an input the next time a simulation is run.**
The active flag is stored in the HDF5 file.

Similarly, results groups are instances of :class:`.ResultNode` class and one of them is always flagged active.
Whenever a simulation is ran, the newly created group is always set as active.

In addition to being used in the next simulation, the active inputs are also used in post-processing and in GUI, so pay attention to what groups are active at the moment.

A group can be set active with

.. code-block:: python

   a5.data.bfield.active.activate()

Each group has a QID which is a 10-number identifier to separate that group from others and it is automatically created when the group data is written to the HDF5 file.
The QID can be used as a reference to the group.

.. code-block:: python

   a5.data.bfield.active.get_qid() # Returns a string with a format "qXXXXXXXXXX"
   a5.data.bfield.q0123456789      # DataGroup corresponding to magnetic field input that has QID = 0123456789

Each :class:`.DataGroup` also has a date when it was created and an user-given description.
The first word (or ten first characters of that word) of the description can be used as a **tag** to refer to the group.

.. code-block:: python

   a5.data.get_date() # Returns date when this group (run) was created.
   a5.data.set_desc("Alpha particle simulation")
   a5.data.ALPHA # Refers to the group on the previous line.

Finally, the treeview can be accessed both attribute-like and dictionary-like.

.. code-block:: python

   # These are all equivalent references to a group
   a5.data.bfield.q0123456789
   a5.data["bfield"]["q0123456789"]
   a5.data.bfield["TAG"]

See the tutorial and API for the following entries for more details.

.. autosummary::
   :nosignatures:

   ~a5py.ascot5io.Ascot5IO
   ~a5py.ascot5io.RunGroup
   ~a5py.ascot5io.InputGroup
   ~a5py.ascot5io.coreio.treedata.DataGroup

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

Neutral beam injectors ``nbi``
==============================

NBI input is used by BBNBI5 and not by ASCOT5.
Hence, its presence is not required when running ASCOT5 simulations.

.. autosummary::
   :nosignatures:

   nbi
   nbi.NBI
   nbi.NBI.write_hdf5

Markers ``marker``
==================

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

===============
Post-processing
===============

Ascotpy
=======

.. currentmodule:: a5py

Many post-processing (and also some pre-processing) tools make use of the Python interface to ASCOT5, referred as ``Ascotpy``.
Instances of :class:`Ascot` automatically establish the interface if ``libascot.so`` has been compiled and found in ``LD_LIBRARY_PATH``.
The interface provides direct access to the C routines and it can be used to interpolate inputs exactly as they are interpolated during the simulation.
It can also be used to simulate markers directly from Python.

The interface is used by first initializing the data (note that this can be memory intensive with large inputs).
The :meth:`~Ascot.input_init` method reads the data from the HDF5 file and initializes it using :mod:`ctypes` so that the data can be passed on to the C routines.
The resources are freed with :meth:`~Ascot.input_free`.
An exception is raised if trying to evaluate a quantity that needs input which has not been initialised.

Here's a complete example on how to initialize input data, evaluate and plot it, and free the resources.

.. code-block:: python

   a5 = Ascot("ascot.h5")
   rgrid = np.linspace(4,8,100)
   zgrid = np.linspace(-2, 2, 100)
   phi   = 0
   t     = 0

   # Initialize currently active bfield
   a5.input_init(bfield=True)

   # Evaluate rho and psi at given (R,phi,z,t) locations
   rho, psi = a5.input_eval(rgrid[0], phi, zgrid[0], t, "rho", "psi")

   # Evaluate rho at given (R,phi,z,t) grid (returns 4D matrix)
   psirz = a5.input_eval(rgrid, phi, zgrid, t, "rho", grid=True)

   # Plot the data (on log10 scale)
   a5.input_plotrz(rgrid, zgrid, "log psi", phi=phi, t=t)

   # Free all resources
   a5.input_free()

   # Initialize bfield with a given QID
   a5.input_init(bfield="0123456789")

   # Free just the bfield data.
   a5.input_free(bfield)

   # Initialize all inputs used by the active run
   a5.input_init(run=True)
   a5.free()

   # Initialize magnetic field input of a run with given QID
   a5.input_init(run="0123456789", bfield=True)
   a5.free()

Complete list of quantities that can be evaluated can be found with :meth:`~Ascot.input_eval_list`.

Methods relevant for input evaluation are listed below.

.. autosummary::
   :nosignatures:

   Ascot
   Ascot.input_init
   Ascot.input_free
   Ascot.input_eval
   Ascot.input_plotrz
   Ascot.get_plasmaquantities
   Ascot.input_eval_list
   Ascot.input_eval_collcoefs
   Ascot.input_getplasmaspecies
   Ascot.input_rhotheta2rz

Simulation output
=================

.. currentmodule:: a5py.ascot5io

The simulation output depends on what diagnostics where active.
Marker ini- and endstate are always present and they contain the marker phase-space position before the first time-step and after the last time-step, respectively.

.. note::
   Plotting routines accept arguments for :obj:`~matplotlib.axes.Axes` where the plot is created (and  axes for the colorbar if applicable).
   This can be used to create figures for publishing as the axes can be customized and reused.
   The axes argument is optional; a new figure is created if the argument is not supplied.

Ini- and endstate
*****************

The data in marker ini- and endstate can be accessed with :meth:`~RunGroup.getstate` method which returns the queried quantity as an array sorted by marker ID.
The marker state contains both particle and guiding-center phase space coordinates, and it also automatically evaluates some derived quantities.
Note that evaluating some of the quantities requires input data to be initialized.
Whether a given quantity is evaluated in particle or guiding-center phase space is dictated by ``mode`` parameter.
Name of the queried quantity must be preceded by "ini" or "end" to indicate whether the quantity is evaluated at the ini- or endstate.
Multiple quantities can be evaluated simultaneously and it is more efficient to do so than to call :meth:`~RunGroup.getstate` separately for each quantity.
Full list of available quantities is given by :meth:`~RunGroup.getstate_list`.

.. code-block:: python

   # Final (R,z) coordinates of marker guiding center and particle positions
   a5.data.active.getstate("end r", "end z", mode="gc")
   a5.data.active.getstate("end r", "end z", mode="prt")

   # Initial guiding center psi coordinate (requires magnetic field to be initialized)
   a5.data.active.getstate("ini psi")

   # Markers can be pruned by their ID and end condition, and the returned quantity
   # can be a difference/relative difference between ini- and ednstate.
   a5.data.active.getstate("diff ekin", endcond="TLIM", ids=[1, 2])

Marker states can be visualized conveniently with :meth:`~RunGroup.plotstate_scatter` and :meth:`~RunGroup.plotstate_histogram`.

.. code-block:: python

   # Quering quantities acts in a similar way as in getstate()
   a5.data.active.plotstate_scatter("ini x", "end y", "diff z", cc="log reldiff ekin",
                                    xmode="prt", endcond="TLIM", ids=[1, 2])

   # 2D scatter plot with a single color
   a5.data.active.plotstate_scatter("end x", "end y", cc="black")

   # 2D scatter plot with color as a coordinate
   a5.data.active.plotstate_scatter("end x", "end y", cc="end ekin")

   # 3D scatter plot with color as a coordinate
   a5.data.active.plotstate_scatter("end x", "end y", "end z", cc="end ekin")

   # 1D stacked histogram where markers with different end conditions are separated
   a5.data.active.plotstate_histogram("ini rho")

   # 2D histogram
   a5.data.active.plotstate_histogram("end rho", "end phimod")

.. autosummary::
   :nosignatures:

   RunGroup
   RunGroup.getstate
   RunGroup.getstate_list
   RunGroup.getstate_markersummary
   RunGroup.plotstate_scatter
   RunGroup.plotstate_histogram


Marker orbits
*************

The orbit output consists of marker phase-space coordinates collected on several points along the marker trajectory.
The data is either gathered at given intervals or it is collected when specific toroidal, poloidal or radial surfaces are crossed.
The latter case allows one to plot Poincaré plots from the orbit data.

The orbit output differs from marker state in that only the phase-space coordinates corresponding to the simulation mode are stored, i.e., GO simulations store only particle coordinates and GC simulations only the guiding-center coordinates.
Other than that, the orbit output is accessed and plotted in same manner as the state output.

.. code-block:: python

   # Return coordinates that are first sorted by marker ID and then by mileage.
   rho, dekin = a5.data.active.getorbit("rho", "diff ekin", endcond=None, ids=None)

   # Plot the trajectory in 2D+1 and 3D+1 where color is the extra dimension
   a5.data.active.plotorbit_trajectory()

.. autosummary::
   :nosignatures:

   RunGroup.getorbit
   RunGroup.getorbit_list
   RunGroup.plotorbit_trajectory
   RunGroup.plotorbit_poincare

End conditions
**************

Marker simulation ends when one or more of the active end conditions is met.
The marker end condition is stored as a bitarray, where the active bits indicate what were the end conditions resulting for that marker's termination.
The post-processing routines use more human-readable approach where end conditions are expressed with strings.
Routines, that allows selecting markers according to their end condition, takes ``endcond`` string as an argument.
This string can either be a single end condition, e.g. ``endcond="WALL"``, in which case all markers **with at least that end condition** are selected.
Providing multiple end conditions, e.g. ``endcond="POLMAX TORMAX"``, selects all markers that have at least those end conditions.
Finally, selecting markers that don't have a given end condition can be done with ``NOT`` keyword as ``endcond="NOT WALL"``.

The following list contains all possible end conditions.

.. autosummary::
   :nosignatures:

   ̃a5py.ascot5io.state.State.ABORTED
   ̃a5py.ascot5io.state.State.NONE
   ̃a5py.ascot5io.state.State.TLIM
   ̃a5py.ascot5io.state.State.EMIN
   ̃a5py.ascot5io.state.State.THERM
   ̃a5py.ascot5io.state.State.WALL
   ̃a5py.ascot5io.state.State.RHOMIN
   ̃a5py.ascot5io.state.State.RHOMAX
   ̃a5py.ascot5io.state.State.POLMAX
   ̃a5py.ascot5io.state.State.TORMAX
   ̃a5py.ascot5io.state.State.CPUMAX
   ̃a5py.ascot5io.state.State.NEUTR
   ̃a5py.ascot5io.state.State.IONIZ

If marker was aborted, something went wrong in the simulation.
**The cause should always be investigated** as the underlying error could have also affected the markers that finished normally.
In most cases the reason for the error is that the marker input contains non-physical data, the input data grids don't cover the whole plasma or the magnetic field ``rho`` is not properly normalized.

To list all end conditions and possible errors, use the :meth:`~RunGroup.getstate_summary` method.
The error shows the error message and the line and file in which it originated.

Distributions
=============

ASCOT5 can be used to collect N-dimensional histograms that represent distribution of the simulated particle population in phase-space.
There are in total of five distributions that can be collected:

  - ``5d`` Distribution in (R, phi, z, ppar, pperp).
  - ``6d`` Distribution in (R, phi, z, pR, pphi, pz).
  - ``rho5d`` Distribution in (rho, theta, phi, ppar, pperp).
  - ``rho6d`` Distribution in (rho, theta, pR, phi, pz).
  - ``com`` Distribution in (mu, ekin, ptor), i.e. in constants of motion space (first adiabatic invariant, energy, canonical toroidal angluar momentum).

While all distributions can be collected irrespective of the simulation mode, there are few considerations.
If you are only interested in the distribution function itself, use either ``5d``/``rho5d`` or ``6d``/``rho6d`` depending on whether the simulation is in guiding -center (5D) or gyro-orbit (6D) mode and in what real-space basis you need the distribution in.
``com`` is for special cases.
However, if you are also interested in moments of the distribution, choose one of the 5D distributions: ``rho5d`` can be used to calculate 1D profiles of the moments while ``5d`` can be used for (R, z) profiles.

In addition to the phase-space coordinates stated here, all distributions (except ``com``) can also have charge state and time as additional coordinates.
Note though that other than that there is no way to separate contributions from different particles.
It is advisable to run separate simulations when there are e.g. several NBIs so that their contributions can be separated.
Finally, all markers should have same mass when distributions are collected as the mass is not a separate coordinate.

Distribution data is obtained in post-processing with the method :meth:`~RunGroup.getdist`, which returns a :class:`~a5py.ascot5io.dist.DistData` object.
For 5D distributions it is possible to change the momentum space from parallel and perpendicular momentum to energy and pitch.

.. code-block:: python

   # Return 5d distribution
   dist = a5.data.active.getdist("5d")

   # 5D distribution where momentum space is converted to energy-pitch
   exidist = a5.data.active.getdist("5d", exi=True)

   # List containing names of the abscissae
   print(dist.abscissae)

   # Get abscissa and edges of an abscissa
   dist.abscissa("r")
   dist.abscissa_edges("r")

   # Distribution (function) and histogram (particles per bin)
   dist.distribution()
   dist.histogram()

The resulting distribution data can be easily sliced, integrated and interpolated using the methods of :class:`~a5py.ascot5io.dist.DistData`.
One can use these methods to process the distribution until it has only one or two non-singleton dimensions left and it can be visualized.

.. code-block:: python

   dist = a5.data.active.getdist("5d")

   # Integrate over time dimension. Also integrate over the charge dimension
   # but use charge as a weight.
   dist.integrate(time=np.s_[:], charge=dist.abscissa("charge"))

   # Take slices and interpolate. These operations won't reduce dimensionality but
   # they reduce the size of the dimensions to one.
   dist.slice(ppar=np.s_[0], pperp=np.s_[0])
   dist.interpolate(phi=2*unyt.deg)

   # Two non-singleton dimensions (R and z) remain so the distribution can be visualized
   a5.data.active.plotdist(dist)

Moments are calculated with :meth:`~RunGroup.getdist_moments` and they are represented with :class:`~a5py.ascot5io.dist.DistMoment` object.
Note that some of the moments require that the input data such as the magnetic field is initialized first.

.. code-block:: python

   # Calculate and plot toroidally averaged density on (R,z) plane
   dist = a5.data.active.getdist("5d")
   mom = a5.data.active.getdist_moments(dist, "density", "chargedensity")
   a5.plotdist_moments(mom, "density")

   # Calculate and plot 1D profile of a toroidally and poloidally averaged density
   dist = a5.data.active.getdist("rho5d")
   mom = a5.data.active.getdist_moments(dist, "density", "chargedensity")
   a5.plotdist_moments(mom, "density")

.. autosummary::
   :nosignatures:

   RunGroup.getdist
   RunGroup.getdist_moments
   RunGroup.plotdist
   RunGroup.plotdist_moments

.. currentmodule:: a5py.ascot5io.dist

.. autosummary::
   :nosignatures:

   DistData
   DistData.distribution
   DistData.histogram
   DistData.abscissa
   DistData.abscissa_edges
   DistData.integrate
   DistData.slice
   DistData.interpolate

   DistMoment
   DistMoment.ordinate

Exceptions
==========

.. currentmodule:: a5py.ascot5io

These are the custom exceptions and warnings used by :mod:`a5py`.

.. autosummary::
   :nosignatures:

   ~a5py.exceptions.AscotIOException
   ~a5py.exceptions.AscotNoDataException
   ~a5py.exceptions.AscotInitException
   ~a5py.exceptions.AscotUnitWarning

Interactive simulations
=======================

Simulations can also be run directly from Python via the interface provided by ``libascot.so`` and :mod:`Ascotpy`.
These simulations are equivalent to running ``ascot5_main`` with the exception that the output is not stored in ``ascot.h5``.
These *interactive* simulations are intended for simulating a handful of markers at most which is why the distribution output is not implemented.
The interface is useful for visualizing marker trajectories or Poincaré plots and evaluating quantities that require markers to be traced only for few orbits.

Before an interactive simulation can be launched, one has to initialize and pack inputs, and set options and markers.

.. code-block:: python

   a5 = Ascot("ascot.h5")

   # You may use any options and markers but in this example we take them from
   # the input file.
   mrk = a5.data.marker.active.read()
   opt = a5.data.options.active.read()

   a5.simulation_initinputs()
   a5.simulation_initmarkers(**mrk)
   a5.simulation_initoptions(**opt)

The marker and options data are dictionaries with same format as required by their corresponding ``write_hdf5`` routines.
The :meth:`Ascot.simulation_initinputs` method requires that no inputs are initialized before calling it.
It differs from the :meth:`Ascot.input_init` method in that it initializes all inputs required for the simulation (except options and markers) and *packs* them.
Packing means that input data is allocated in a single array and this is what ASCOT5 does internally when the data is *offloaded* for simulation.
When inputs are packed, the input data cannot be deallocated or altered but otherwise it is possible to use methods like :meth:`Ascot.input_eval`, that evaluate and plot input data, normally.

Simulation is run and data is accessed as:

.. code-block:: python

   # Running simulations returns a virtual run that acts the same way as normal runs
   vrun = a5.simulation_run()
   r, z = vrun.getorbit("r", "z")

   # Simulations can be rerun by freeing the previous results first
   # Options can be changed when results are freed
   a5.simulation_free(diagnostics=True)
   a5.simulation_initoptions(**opt)
   vrun = a5.simulation_run()

   # Simulations can be rerun with new markers if markers are freed as well
   a5.simulation_free(markers=True, diagnostics=True)
   a5.simulation_initmarkers(**mrk)
   vrun = a5.simulation_run()

   # Free all resources
   a5.simulation_free()

List of methods relevant for running live simulations can be found below.

.. autosummary::
   :nosignatures:

   a5py.Ascot.simulation_initinputs
   a5py.Ascot.simulation_initmarkers
   a5py.Ascot.simulation_initoptions
   a5py.Ascot.simulation_run
   a5py.Ascot.simulation_free

Units
=====

:mod:`a5py` uses :mod:`unyt` package to include units in physical quantities.

Unyt package is easy to use as the documentation shows and briefly summarized here.

.. code-block:: python

   # Units can be assigned to values or numpy arrays by multiplying with the unit
   a = 1 * unyt.m
   b = np.array([1, 2]) * unyt.s

   # Units obey basic arithmetic
   velocity = a / b
   velocity.units # Has units m/s

   # Strip units
   velocity.v

   # Assign units using a string
   energy = unyt.unyt_quantity.from_string("10 eV")

   # Converting to different units (returns a new array)
   energy_J = energy.to("J")

   # Converting in place
   energy_J.convert_to_units("eV")

   # Converting to SI units
   energy.convert_to_base("mks")

   Converting to "ascot" unit system
   energy.convert_to_base("ascot")

The package takes care of unit conversions automatically when operating quantities that have different units.
Exception is raised whenever trying to operate quantities that don't have matching physical dimensions.

:mod:`a5py` implements its own unit system ``ascot`` which for the most part is in SI units with the exception of energy [eV], mass [amu], charge [e], and angle [deg].

Functions expecting arguments to have units issue warning if units are not provided and assumes that the arguments already are in correct units.
The warning can be suppressed with:

.. code-block:: python

   import warnings
   from a5py import AscotUnitWarning
   warnings.filterwarnings("ignore", category=AscotUnitWarning)
