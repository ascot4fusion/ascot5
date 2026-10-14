.. _Data:

====
Data
====

.. admonition:: Summary

   - All ASCOT5 inputs and results are contained in a single HDF5 file.
   - Data in the HDF5 should not be directly accessed and the Python API should be used instead.
   - Each input and each result has a ten number quasi-unique identifier (QID) randomly assigned when the data is written.
   - Each result group stores QIDs of the inputs that were used in the simulation.

The ASCOT5 HDF5 file, named ``ascot.h5`` from here on out, contains all data required to run a simulation.
All data is contained in a single file so that the results, and the inputs that produced them, are always kept together.
The file is designed to hold multiple inputs (even of same type) and results of multiple simulations.
One is encouraged to use a single file for the whole project or study.

The exact structure of the contents of the HDF5 file is not that relevent since it should always be accessed via the Python interface provided by :mod:`a5py`.
When accessed via Python, the contents of the file as they appear are illustrated below:


ASCOT object
============


Options
=======


Command line tools
==================

Some useful scripts can be found in `./bin` folder or accessed directly from terminal if ``a5py`` was installed succesfully.
These scripts modify or disbplay the contents of the HDF5 file and they are listed here.

.. list-table:: Command line scripts
   :widths: 25 75
   :header-rows: 0

   * - ``a5manage``
     - Manage data groups by setting them as active, removing them, or copying them to another file.
       Inputs that have been used by a run cannot be removed unless that run is removed first.
   * - ``a5gui``
     - Open GUI where you can manage your data and view results and inputs.

End conditions
==============

AFSI
====

AFSI calculates distribution of fusion products from the interaction of two arbitrary populations.
It can be used to create either a fast particle source for ASCOT5 simulations or a neutron source e.g. for codes such as Serpent.
AFSI has three modes of operation: *thermal* where two Maxwellian populations are interacting, *beam-thermal* where one population is Maxwellian and the other is given by arbitrary distribution (e.g. beam ion slowing-down distribution), or *beam-beam* where both populations are given by arbitrary distributions.
The distributions must have the format of a 5D distribution.

AFSI uses ``libascot.so`` to perform the calculations efficiently, but the interface is completely in Python.
Each AFSI5 run is stored in HDF5 and its Python counterpart, , operates similarly to .
The main difference is that the group contains only distribution data: ``prod1`` and ``prod2`` which are :class:`~a5py.ascot5io.dist.DistData` objects representing the fusion products.

Running AFSI requires that the inputs contain magnetic field and plasma data.
AFSI can automatically detect the reactants from the plasma input when running thermal or beam-thermal simulations.

.. currentmodule:: a5py.routines.afsi5

To run AFSI, use :class:`.Ascot` attribute ``afsi`` which is an instance of :class:`Afsi`.
Here's a demonstration of running AFSI withh different modes:

.. code-block:: python

   a5 = Ascot("ascot.h5")
   # Lists all available reactions
   a5.afsi.reactions()

   # Thermal mode requires defining the abscissae for the product distributions.
   # For the other modes the output distributions have the same size and dimensions
   # as the beam input distribution.
   a5.afsi.thermal(
      "DT_He4n",
      rmin, rmax, nr, zmin, zmax, nz,
      minphi=0, maxphi=2*np.pi, nphi=1, nmc=1000,
      minppara=-1.3e-19, maxppara=1.3e-19, nppara=80,
      minpperp=0, maxpperp=1.3e-19, npperp=40)

   # In beam-thermal, the beam is the first reactant unless swap=True
   a5.afsi.beamthermal("DD_Tp", beamdist, swap=False)

   # Beam-beam can either be two separate beams or the beam interacting with itself
   a5.afsi.beambeam("DD_Tp", beamdist, beam2=None)

   # Results are accessed as usual
   dist = a5.data.active.getdist("prod1")

.. autosummary::
   :nosignatures:

   Afsi
   Afsi.thermal
   Afsi.beamthermal
   Afsi.beambeam
   Afsi.reactions

BBNBI
=====

BBNBI implements a beamlet based model for generating neutrals from an injector geometry.
The code follows generated neutrals until they are either ionized or they hit the wall (and become shinethrough).

BBNBI is a separate program that has to be compiled separately: ``make bbnbi5`` in the main folder.
This creates binary ``./build/bbnbi5`` that uses the same input file and shares some inputs with the main program.
A ``bbnbi5`` run requires following inputs: ``nbi``, ``bfield``, ``plasma``, ``wall``, and ``options``.
Options are only used partially as only the distribution options are used to bin ionized markers to create particle source distribution.
The ``nbi`` input refers to :class:`.NBI` which consists of a single or a bundle of injectors represented by :class:`.Injector`.

.. note::
   It is important to note that there are now way to separate results by injector if the injectors were bundled.
   Only bundle injectors if they have identical geometries or otherwise there is no reason to separate the results.

For testing purposes one can use to create a simple injector and there are scripts to generate injectors for various machines, so generating the injector data explicitly is rarely necessary.
To generate injectors explicitly, refer to :class:`.Injector` for the required data.
In short, one has to provide the location of the beamlets, their direction and beam divergence, injected species, the energy fractions and the beam power.

.. note::
   Scripts to generate injector geometries exist for ITER, JET, ASDEX Upgrade, MAST-U, JT-60SA, NSTX, and TCV.
   These are stored in a private repository so contact the ASCOT team in order to gain access.

Once all the required inputs are present, BBNBI is run in a same way as ASCOT5 with the exception that BBNBI is not MPI parallelized:

.. code-block:: bash

   ./bbnbi5 --in=ascot.h5 --d="MYTAG My description" --n=100000 --t0=0.0 --t1=1.0

Here the command line arguments specify the number of markers to be injected, ``n``, and the time-interval ``[t0,t1]`` when the beam is enabled.
The time interval is only relevant for time-dependent simulations and all it does is to randomly pick a birth time for an injected marker from the interval.

Each BBNBI run is stored in HDF5 and the data is accessed via :class:`.BBNBIGroup` which in many ways is similar to :class:`.RunGroup`.
Instead of storing separate ini- and endstates, only the final state of the markers is stored and it can be used to evaluate wall loads in post-processing or to initialize markers for the follow-up ASCOT5 simulation.

.. code-block:: python

   # File which has BBNBI5 run tagged "BBNBI"
   a5 = Ascot("ascot.h5")
   # List all inputs and outputs
   a5.data.BBNBI.ls()
   # Marker summary and data (note that there is only one "state")
   a5.data.BBNBI.getstate_markersummary()
   r, z = a5.data.BBNBI.getstate("r", "z")
   # Shinethrough
   a5.input_init(bfield=True)
   a5.data.active.plotwall_torpol()
   a5.input_free()
   # Create marker input from ionized markers
   ids = a5.data.active.getstate("ids", endcond="IONIZED")
   mrk = a5.data.active.getstate_markers("gc", ids=ids)
   a5.data.create_input("gc", **mrk)
   # Read the particle source distribution
   dist = a5.data.getdist("5d")

.. currentmodule:: a5py.ascot5io

.. autosummary::
   :nosignatures:

   BBNBIGroup
   BBNBIGroup.getstate
   BBNBIGroup.getstate_markers
   BBNBIGroup.getstate_markersummary
   BBNBIGroup.getdist
   BBNBIGroup.getwall_figuresofmerit
   BBNBIGroup.getwall_loads
   BBNBIGroup.getwall_3dmesh
   BBNBIGroup.plotwall_3dstill
   BBNBIGroup.plotwall_3dinteractive
   BBNBIGroup.plotwall_torpol

BioSaw
======

BioSaw is a solver using Biot-Savart law to conveniently include perturbation from external coils to ASCOT5 simulation.
It computes the magnetic field from a coil geometry and turns it into a :class:`.B_3DS` input.
The produced data is divergence-free although the level of divergence in simulation, when the data is interpolated with splines, depends on the grid resolution.
BioSaw uses ``libascot.so`` to perform the calculations efficiently, but the interface is completely in Python.

To run BioSaw, you need a file that has either :class:`.B_2DS` or :class:`.B_3DS` magnetic field present or set as active.
One starts by defining the coil geometry:

.. code-block:: python

   rad   = 3.0
   ang   = np.linspace(0, 2*np.pi, 40)
   coilx = rad * np.cos(ang) + 6.2
   coily = np.zeros(ang.shape)
   coilz = rad * np.sin(ang)
   coilxyz = np.array([coilx, coily, coilz])

.. note::
   The coil geometry does not have to form a closed loop as BioSaw simply calculates the field arising from each line segment.

.. currentmodule:: a5py.routines.biosaw5

BioSaw is represented by class :class:`~a5py.routines.biosaw5.BioSaw` whose instance exists as an attribute named ``biosaw`` in :class:`.Ascot`.
The method :meth:`BioSaw.calculate` can be used to calculate the magnetic field explicitly, but the methods :meth:`BioSaw.addto2d` and :meth:`BioSaw.addto3d` are more convenient to use when one wants to include the perturbation to an existing 2D or 3D field.

To include perturbation from a coil to field that is already 3D:

.. code-block:: python

   # Current is in Ampéres
   b3d = a5.biosaw.addto3d(coilxyz.T, current=30e4, revolve=None)

   # Multiple coils can be given as lists
   b3d = a5.biosaw.addto3d([coilxyz.T, ...], current=[30e4, ...], revolve=None)

In 2D, the most common application is the inclusion of TF ripple (which overrides the 2D toroidal field already present in the data).
Here it is not necessary to know the coil current, as the current in TF coils can be scaled so that the average toroidal magnetic field on axis has the given value.
Furthermore when multiple identical coils has a toroidally symmetrical arrangement, the field can be calculated with just a single coil and the total field computed by revolving and summing this result.
This is achieved with ``revolve`` which defines the number of toroidal grid points (not the toroidal angle) the field is revolved in each rotation until we have gone whole turn.
So if there are e.g. 18 TF coils, the number of grid points must be divisible by 18 and revolve is the number of grid points divided by the number of coils.

.. code-block:: python

   # Since the original field is 2D, one can now define the toroidal resolution for the 3D field
   b3d = a5.biosaw.addto2d(
         coilxyz.T, phimin=0, phimax=360, nphi=180,
         revolve=10, b0=5.3)

   # Write the data to HDF5
   a5.data.create_input("B_3DS", **b3d, desc="TFRIPPLE")

.. autosummary::
   :nosignatures:

   BioSaw
   BioSaw.addto2d
   BioSaw.addto3d
   BioSaw.calculate


Data access
===========


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


.. code-block:: python

   import warnings
   from a5py import AscotUnitWarning
   warnings.filterwarnings("ignore", category=AscotUnitWarning)

Note on units
=============

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

Functions expecting arguments to have units issue warning if units are not provided and assumes that the arguments already are in correct units.
The warning can be suppressed with:
