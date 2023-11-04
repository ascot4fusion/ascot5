.. _Postprocessing:

===============
Post-processing
===============

.. currentmodule:: a5py

.. admonition:: Summary

   - While some postprocessing can be done without ``ascotpy``, the Python interface to ASCOT5, most advanced features require it.

     - ``ascotpy`` is enabled when the C library ``libascot.so`` has been compiled.

   - Input data can be interpolated and visualized using the ``input_*`` methods of :class:`Ascot`.
   - The run group, :class:`~ascot5io.RunGroup`, has methods for accessing and visualizing the output.

     - ``a5.data.active.getstate(*quantities, mode=..., ids=..., endcond=...)`` fetch marker ini/endstate.
     - ``a5.data.active.getorbit(*quantities, ids=..., endcond=...)`` fetch orbit data.
     - ``a5.data.active.getdist(dist)`` fetch distribution data and ``a5.data.active.getmoments(dist)`` calculates the moments.
     - ``a5.data.active.getwall`` evaluates wall loads.

   - What output is present depends on what diagnostics were enabled in the simulation options.
   - Simulations can be run directly from Python using the ``simulate_*`` methods of :class:`Ascot`.
   - Supporting tools BioSaw and AFSI are run via ``a5.biosaw`` and ``a5.afsi`` objects, respectively, while BBNBI is a separate program that must be compiled and run similarly as ASCOT5.

Simulation output is accessed via :class:`~ascot5io.RunGroup` object, i.e. `a5.data.active` is by default the object containing the recent simulation results.
Its contents depend on what diagnostics were enabled in the simulation options, except marker ini- and endstate (marker phase space position and associated quantities at the beginning and at the end of the simulation) which are always present.

Many post-processing (and also some pre-processing) tools make use of the Python interface to ASCOT5, referred as ``ascotpy``.
Using ``ascotpy`` requires that ``libascot.so`` has been compiled as per the installation instructions.

.. note::
   Most plotting routines accept arguments for :obj:`~matplotlib.axes.Axes` where the plot is created (and  axes for the colorbar if applicable).
   This can be used to create figures for publishing as the axes can be customized and reused.
   The axes argument is optional; a new figure is created if the argument is not supplied.

Using the Python interface
==========================

Instances of :class:`Ascot` automatically establish the interface if ``libascot.so`` has been compiled.
The interface provides direct access to the C routines and it can be used to interpolate inputs exactly as they are interpolated during the simulation.
It can also be used to simulate markers directly from Python.

To use the interface, start by initializing the input data with :meth:`~Ascot.input_init`:

.. code-block:: python

   a5 = Ascot("ascot.h5")

   # Initialize currently active bfield
   a5.input_init(bfield=True)

.. note::
   If the ``VERBOSE`` level is not zero, ASCOT5 prints information on the initialized inputs to stdout.
   Same is true for ``libascot.so`` but in the Python interface these messages are supressed.
   These messages can be made visible again with ``a5 = Ascot("ascot.h5", mute="no")`` which can be helpful for investigating when there are problems with inputs.

This function reads the input data from the HDF5 file, and passes it to C-routines where it is initialized and stored in the memory.
After initialization, one can use :meth:`~Ascot.input_eval` to interpolate the data using the very same C-routines which ASCOT5 uses when it traces markers:

.. code-block:: python

   import unyt # units are not necessary but otherwise warnings are issued
   a5 = Ascot("ascot.h5")
   rgrid = np.linspace(4,8,100)*unyt.m
   zgrid = np.linspace(-2, 2, 100)*unyt.m
   phi   = 0*unyt.deg
   t     = 0*unyt.s

   # Evaluate rho and psi at given (R,phi,z,t) locations
   rho, psi = a5.input_eval(rgrid[0], phi, zgrid[0], t, "rho", "psi")

   # Evaluate rho at given (R,phi,z,t) grid (returns 4D matrix)
   psirz = a5.input_eval(rgrid, phi, zgrid, t, "rho", grid=True)

Complete list of quantities that can be evaluated can be found with :meth:`~Ascot.input_eval_list`.
Note however, that :meth:`~Ascot.input_eval` is just one of the many functions that use ``ascotpy``.
Once you are done, free the resources with :meth:`~Ascot.input_free`:

.. code-block:: python

   a5.input_free()

An exception is raised if one tries to call a routine that expects specific inputs to be initialized.
It is in user's control to decide what input is used and it doesn't have to be the same as was used in the simulation.
Below is an example on different ways to specify what input is initialized:

.. code-block:: python

   # Init active input
   a5.input_init(bfield=True)

   # Init input with a given QID
   qid = a5.data.bfield.TAG.get_qid()
   a5.input_init(bfield=qid)

   # Init all inputs used by a given run
   a5.input_init(run=True)

   # Init just the bfield of a run with a given QID
   qid = a5.data.TAG.get_qid()
   a5.input_init(run=qid, bfield=True)

.. note::
   Initializing inputs, especially the magnetic field, can be memory-intensive since this is usually done on a desktop/laptop and not in a dedicated computing node.
   However, while simulations require magnetic field with fine resolution, using one with a coarser grid, and thus less memory-intensive, for post-processing is acceptable.

.. rubric:: Classes and methods for input processing

.. autosummary::
   :nosignatures:

   Ascot
   Ascot.input_init
   Ascot.input_free
   Ascot.input_eval
   Ascot.input_plotrz
   Ascot.input_eval_list
   Ascot.input_rhotheta2rz
   Ascot.input_rhovolume
   Ascot.input_eval_safetyfactor
   Ascot.input_eval_ripple
   Ascot.input_plotwallcontour
   Ascot.input_getplasmaspecies
   Ascot.input_eval_collcoefs
   Ascot.input_eval_atomicsigma

.. currentmodule:: a5py.ascot5io

End conditions
==============

Marker simulation ends when one or more of the active end conditions is met.
The marker end condition is stored as a bit array, where the active bits indicate what were the end conditions resulting for that marker's termination.
The post-processing routines use more human-readable approach where end conditions are expressed with strings.
This string can either be a single end condition, e.g. ``endcond="WALL"``, in which case all markers **with at least that end condition** are selected.
Providing multiple end conditions, e.g. ``endcond="[POLMAX, TORMAX]"``, selects all markers that have at least all of those end conditions.
Finally, selecting markers that don't have a given end condition can be done with ``NOT`` keyword as ``endcond="NOT WALL"``.

The following list contains all possible end conditions.

.. autosummary::
   :nosignatures:

   ~a5py.ascot5io.state.State.ABORTED
   ~a5py.ascot5io.state.State.NONE
   ~a5py.ascot5io.state.State.TLIM
   ~a5py.ascot5io.state.State.EMIN
   ~a5py.ascot5io.state.State.THERMAL
   ~a5py.ascot5io.state.State.WALL
   ~a5py.ascot5io.state.State.RHOMIN
   ~a5py.ascot5io.state.State.RHOMAX
   ~a5py.ascot5io.state.State.POLMAX
   ~a5py.ascot5io.state.State.TORMAX
   ~a5py.ascot5io.state.State.CPUMAX
   ~a5py.ascot5io.state.State.HYBRID
   ~a5py.ascot5io.state.State.NEUTRAL
   ~a5py.ascot5io.state.State.IONIZED

.. warning::
   If marker was aborted, something went wrong in the simulation.
   **The cause should always be investigated** as the underlying error could have also affected the markers that finished normally.
   Most common causes of aborted markers is that the marker input contains non-physical data, the input data grids don't cover the whole plasma or that the magnetic field ``rho`` is not properly normalized.

To list all end conditions present in the simulation output, and possible errors, use :meth:`~RunGroup.getstate_summary`.
The listed errors shows the error message and the line and file in which it originated.

Ini- and endstate
=================

The data in marker ini- and endstate can be accessed with :meth:`~RunGroup.getstate` method which returns the queried quantity as an array sorted by marker ID.
The marker state contains both particle and guiding-center phase space coordinates, and the method also automatically evaluates some derived quantities.
Whether a given quantity is evaluated in particle or guiding-center phase space is dictated by ``mode`` parameter.

Name of the queried quantity must be preceded by "ini" or "end" to indicate whether the quantity is evaluated at the ini- or endstate.
Multiple quantities can be evaluated simultaneously and it is more efficient to do so than to call :meth:`~RunGroup.getstate` separately for each quantity.

Full list of available quantities is given by :meth:`~RunGroup.getstate_list`.
Note that evaluating some of the quantities requires input data to be initialized.

Here are examples on how to query marker states:

.. code-block:: python

   # Final (R,z) coordinates of marker guiding center and particle positions
   a5.data.active.getstate("end r", "end z", mode="gc")
   a5.data.active.getstate("end r", "end z", mode="prt")

   # Initial guiding center psi coordinate (requires magnetic field to be initialized)
   a5.data.active.getstate("ini psi")

   # Markers can be pruned by their ID and end condition, and the returned quantity
   # can be a difference/relative difference between ini- and ednstate.
   a5.data.active.getstate("diff ekin", endcond="TLIM", ids=[1, 2])

Marker states can be visualized conveniently with :meth:`~RunGroup.plotstate_scatter` and :meth:`~RunGroup.plotstate_histogram`:

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

.. rubric:: Marker evaluation and plotting

.. autosummary::
   :nosignatures:

   RunGroup.getstate
   RunGroup.getstate_list
   RunGroup.getstate_markersummary
   RunGroup.plotstate_scatter
   RunGroup.plotstate_histogram


Marker orbits
=============

The orbit output consists of marker phase-space coordinates collected on several points along the marker trajectory.
The data is either *i*) gathered at given intervals or *ii*) collected when specific toroidal, poloidal or radial surfaces are crossed.
The latter case allows one to plot Poincaré plots from the orbit data.

The orbit output differs from marker state in that only the phase-space coordinates corresponding to the simulation mode are stored, i.e., GO simulations store only particle coordinates and GC simulations only the guiding-center coordinates.
Also the returned values are first sorted by marker ID and then by the mileage.
Other than that, the orbit output is accessed and plotted in same manner as the state output.

.. code-block:: python

   # Return rho value and the change in energy for markers that were not
   # lost to the wall. Orbit data can also be parsed by marker ID
   rho, dekin = a5.data.active.getorbit("rho", "diff ekin", endcond="NOT WALL", ids=None)

   # Plot the trajectories of all markers in 2D+1D (color)
   a5.data.active.plotorbit_trajectory("r", "z", c="pitch")

.. rubric:: Marker evaluation and plotting

.. autosummary::
   :nosignatures:

   RunGroup.getorbit
   RunGroup.getorbit_list
   RunGroup.plotorbit_trajectory
   RunGroup.plotorbit_poincare

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

.. note::
   In addition to the phase-space coordinates stated here, all distributions (except ``com``) can also have charge state and time as additional coordinates.
   Other than that there is no way to separate contributions from different particles.
   Therefore, it is advisable to run separate simulations for separate sources, e.g. when there are several NBIs, so that their contributions can be separated.
   Finally, all markers in a simulation should have same mass when distributions are collected as the mass is not a separate coordinate.

Distributions and their moments are represented with :class:`~a5py.ascot5io.dist.DistData` and :class:`~a5py.ascot5io.dist.DistMoment` objects, respectively.
The distribution object is created when :meth:`~RunGroup.getdist` is called, and then the object can be processed (which has no effect on the data on the disk).
For 5D distributions it is possible to change the momentum space from parallel and perpendicular momentum to energy and pitch.

Example illustrating how distributions can be processed:

.. code-block:: python

   # Return 5d distribution
   dist = a5.data.active.getdist("5d")

   # Same 5D distribution but where momentum space is converted to energy-pitch
   exidist = a5.data.active.getdist("5d", exi=True)

   # Distribution object has both methods and attributes that can be used
   dist.abscissae           # List containing names of the abscissae
   dist.abscissa("r")       # Abscissa bin center values
   dist.abscissa_edges("r") # Abscissa bin edges

   dist.distribution() # The physical distribution function
   dist.histogram()    # Distribution histogram (particles per bin)

The resulting distribution data can be easily sliced, integrated and interpolated using the methods of :class:`~a5py.ascot5io.dist.DistData`:

.. code-block:: python

   dist = a5.data.active.getdist("5d")

   # Integrate over time dimension. Also integrate over the charge dimension
   # but use charge as a weight.
   dist.integrate(time=np.s_[:], charge=dist.abscissa("charge"))

   # Take slices and interpolate. These operations won't reduce dimensionality but
   # they reduce the size of the dimensions to one.
   dist.slice(ppar=np.s_[0], pperp=np.s_[0])
   dist.interpolate(phi=2*unyt.deg)

One can use these methods to process the distribution until it has only one or two non-singleton dimensions left so that it can be visualized:

.. code-block:: python

   dist.plot()

Moments are calculated with :meth:`~RunGroup.getdist_moments`.
Note that some of the moments require that the input data such as the magnetic field is initialized first.

Example on how to calculate and visualize moments:

.. code-block:: python

   # Calculate and plot toroidally averaged density on (R,z) plane
   dist = a5.data.active.getdist("5d")
   mom = a5.data.active.getdist_moments(dist, "density", "chargedensity")
   mom.plot("density")

   # Calculate and plot 1D profile of a toroidally and poloidally averaged density
   dist = a5.data.active.getdist("rho5d")
   mom = a5.data.active.getdist_moments(dist, "density", "chargedensity")
   mom.plot("density")

.. rubric:: Retrieving distributions and moments

.. autosummary::
   :nosignatures:

   RunGroup.getdist
   RunGroup.getdist_moments

.. rubric:: Distribution and moment data

.. currentmodule:: a5py.ascot5io.dist

.. autosummary::
   :nosignatures:

   DistData
   DistData.abscissa
   DistData.abscissa_edges
   DistData.distribution
   DistData.histogram
   DistData.phasespacevolume
   DistData.plot
   DistData.integrate
   DistData.slice
   DistData.interpolate

   DistMoment
   DistMoment.ordinate
   DistMoment.add_ordinates
   DistMoment.list_ordinates
   DistMoment.plot

Interactive simulations
=======================

.. currentmodule:: a5py

Simulations can also be run directly from Python via the interface provided by ``ascotpy``.
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
Packing means that input data is allocated in a single array which is what ASCOT5 does internally when the data is *offloaded* for simulation.
Once the inputs are packed, the input data cannot be deallocated or altered but otherwise it is possible to use methods like :meth:`Ascot.input_eval`, that evaluate and plot input data, normally.

Running a simulation returns a :class:`~a5py.routines.virtualrun.VirtualRun` object that acts similarly as an ordinary run object :class:`~a5py.ascot5io.RunGroup` would.

.. code-block:: python

   # Running simulations returns a virtual run that acts the same way as normal runs
   vrun = a5.simulation_run()
   r, z = vrun.getorbit("r", "z")

Since the results are now kept in the memory and not stored in the disk, rerunning a simulation requires that the previous results are free'd first.
Note that the virtual run object cannot be used once its results have been free'd.

.. code-block:: python

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

.. rubric:: Tools for the interactive simulations

.. autosummary::
   :nosignatures:

   Ascot.simulation_initinputs
   Ascot.simulation_initmarkers
   Ascot.simulation_initoptions
   Ascot.simulation_run
   Ascot.simulation_free
   ~a5py.routines.virtualrun.VirtualRun

Marker generation
=================

The marker generation serves two purposes.
First, the generated population must represent the particle population it models.
Second, the marker generation should aim to gain best possible signal-to-noise ratio so that one can use as few markers as possible to obtain the results.
For example, alpha particles are mostly born at the core but simulating the markers at the core makes little sense if we are only interested in losses which mostly originate at the edge.

For testing purposes one can just manually initialize markers at the desired position.
Here is a helpful tool as it initializes arrays of dummy markers of given particle species:


Another helpful tool for testing and debugging purposes is which converts marker end states to new marker input:


To generate markers for actual simulation, one can use AFSI to generate fusion alphas and BBNBI to generate beam ions.
These tools create 5D distributions which can then be given as an input for the marker generator :class:`` that samples markers from the distribution.
The marker generator takes two distributions as an input: one that represents the physical particle population and from which the weights are assigned, and one probability distribution from which markers are sampled.

.. currentmodule:: a5py.ascot5io

.. autosummary::
   :nosignatures:

   RunGroup.getstate_markers

AFSI
====

AFSI calculates distribution of fusion products from the interaction of two arbitrary populations.
It can be used to create either a fast particle source for ASCOT5 simulations or a neutron source e.g. for codes such as Serpent.
AFSI has three modes of operation: *thermal* where two Maxwellian populations are interacting, *beam-thermal* where one population is Maxwellian and the other is given by arbitrary distribution (e.g. beam ion slowing-down distribution), or *beam-beam* where both populations are given by arbitrary distributions.
The distributions must have the format of a 5D distribution.

AFSI uses ``libascot.so`` to perform the calculations efficiently, but the interface is completely in Python.
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

BBNBI
=====

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

BioSaw
======

BioSaw is a solver using Biot-Savart law to conveniently include perturbation from external coils to ASCOT5 simulation.
It computes the magnetic field from a coil geometry and turns it into a :class:`.B_3DS` input.
The produced data is divergence-free although the level of divergence in simulation, when the data is interpolated with splines, depends on the grid resolution.
BioSaw uses ``libascot.so`` to perform the calculations efficiently, but the interface is completely in Python.

A common use of BioSaw is to incorporate the ripple from TF coils into an axisymmetric field.
Begin by defining the coil geometry

.. note::
   The coil geometry does not have to form a closed loop as BioSaw simply calculates the field arising from each line segment.

To run BioSaw, first open a file that has either :class:`.B_2DS` or :class:`.B_3DS` magnetic field present or set as active.
The calculation is then simply done as follows using the ``biosaw`` attribure which is an instance of :class:`.BioSaw`:

.. code-block:: python

   # Define a simple coil which does not have to form a closed loop


   # Open a file that has B_2DS input "B2D" and B_3DS input "B3D"
   # BioSaw is operated via "biosaw" attribute.
   a5 = Ascot("ascot.h5")

   b3d = a5.biosaw.addto2d(
             coilxyz.T, phimin=0, phimax=360, nphi=180,
             revolve=10, b0=5.3)

.. autosummary::
   :nosignatures:

   a5py.routines.biosaw5.BioSaw
   a5py.routines.biosaw5.BioSaw.addto2d
   a5py.routines.biosaw5.BioSaw.addto3d
   a5py.routines.biosaw5.BioSaw.calculate
