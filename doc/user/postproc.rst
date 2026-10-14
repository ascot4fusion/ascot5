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
