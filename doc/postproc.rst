.. _Postprocessing:

===============
Post-processing
===============

.. currentmodule:: a5py

Many post-processing (and also some pre-processing) tools make use of the Python interface to ASCOT5, referred as ``Ascotpy``.
Instances of :class:`Ascot` automatically establish the interface if ``libascot.so`` has been compiled.
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
*************

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

Advanced
========

Interactive simulations
***********************

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

Neutral beam simulations (BBNBI5)
*********************************

BBNBI5 implements a beamlet based model for generating neutrals from an injector geometry.
Neutrals are then followed until they are ionized or they hit the wall (shinethrough).

BBNBI5 uses :class:`.NBI` input which consists of a single or a bundle of injectors represented by :class:`.Injector`.
It is important to note that there are now way to separate results by injector if the injectors were bundled.
Only bundle injectors if they have identical geometries or otherwise there is no reason to separate the results.

BBNBI5 is a separate program that has to be compiled separately: ``make bbnbi5`` in the main folder.
This creates binary ``./build/bbnbi5`` that uses the same input file and shares some inputs with the main program.
A ``bbnbi5`` run requires following inputs: ``nbi``, ``bfield``, ``plasma``, ``wall``, and ``options``.
The only options that are used from the input are settings for the distributions: if a distribution is toggled on, BBNBI5 gathers ionized markers on that distribution.
When these are present, the code is run with

.. code-block::

   # n is the number of markers to be used in simulation
   # writemarkers creates marker input from ionized markers (on by default)
   # [t0, t1] is the time interval when the beam is on. This affects the time
   # assigned to generated markers and the output distributions are weighted
   # with (t1-t0). Default is t0=t1=0.0 where all markers have t=0.0 and the
   # the distribution is weighted with 1.0 s
   ./bbnbi5 --in=ascot.h5 --d="MYTAG My description" --n=1000 --writemarkers=1 --t0=0.0 --t1=1.0

If the option was set to generate marker input, there is a new instance of ``marker`` present in the file.
Each BBNBI5 run is also stored in HDF5, and in Python it can read and access the data in same fashion as the ASCOT5 runs.
Open Python terminal to access the data

.. code-block:: python

   # File which has BBNBI5 run tagged "BBNBI"
   a5 = Ascot("ascot.h5")
   # List all inputs and outputs
   a5.data.BBNBI.ls()

.. autosummary::
   :nosignatures:

   a5py.Ascot.simulation_initinputs

Generating fusion products (AFSI5)
**********************************

AFSI5 calculates distribution of fusion products from the interaction of two arbitrary populations.
It can be used to create either a fast particle source for ASCOT5 simulations or a neutron source e.g. for codes such as Serpent.
AFSI5 has three modes of operation: *thermal* where two Maxwellian populations are interacting, *beam-thermal* where one population is Maxwellian and the other is given by arbitrary distribution (e.g. beam ion slowing-down distribution), or *beam-beam* where both populations are given by arbitrary distributions.
The distributions must have the format of a 5D distribution.

AFSI5 uses ``libascot.so`` to perform the calculations efficiently, but the interface is completely in Python.
Each AFSI5 run is also stored in HDF5, and in Python it can read and access the data in same fashion as the ASCOT5 runs.
Assuming that we have a file with magnetic field data and DT plasma input present, and also a run tagged "BEAMS" which contains beam slowing-down distribution in 5D, then AFSI5 is run as

.. code-block:: python

   a5 = Ascot("ascot.h5")
   a5.input_init(bfield=True, plasma=True)
   a5.afsi.
   a5.input_free()

   # The output is stored in HDF5 as a results group (and set as active):
   a5.ls()
   a5.data.active.plot_dist("r", "z", "product1")

Computing 3D field from coils (BioSaw)
**************************************

BioSaw is a solver using Biot-Savart law to conveniently include perturbation from external coils to ASCOT5 simulation.
It computes the magnetic field from a coil geometry and turns it into a :class:`.B_3DS` input.
A common use of BioSaw is to incorporate the ripple from TF coils or the error-field from RMP coils into simulations.
The produced data is divergence-free although the level of divergence in simulation, when the data is interpolated with splines, depends on the grid resolution.

BioSaw uses ``libascot.so`` to perform the calculations efficiently, but the interface is completely in Python.
To run BioSaw, first open a file that has either :class:`.B_2DS` or :class:`.B_3DS` magnetic field present or set as active.
The calculation is then simply done as follows using the ``biosaw`` attribure which is an instance of :class:`.BioSaw`:

.. code-block:: python

   # Define a simple coil which does not have to form a closed loop


   # Open a file that has B_2DS input "B2D" and B_3DS input "B3D"
   # BioSaw is operated via "biosaw" attribute.
   a5 = Ascot("ascot.h5")

.. autosummary::
   :nosignatures:

   a5py.routines.biosaw5.BioSaw
   a5py.routines.biosaw5.BioSaw.addto2d
   a5py.routines.biosaw5.BioSaw.addto3d
   a5py.routines.biosaw5.BioSaw.calculate
