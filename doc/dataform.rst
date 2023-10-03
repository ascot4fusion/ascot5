.. _Data:

====
Data
====

Simulation inputs and options, as well as the outputs once the simulation is complete, are stored in a single HDF5 file.

The ASCOT5 HDF5 file, ``ascot.h5`` from here on out, contains all data required to run a simulation, and so the results and the inputs that produced them are always kept together.
The file is designed to hold multiple inputs (even of same type) and results of multiple simulations to encourage one to use a single file for the whole project or study.

The exact structure of the HDF5 file is not that relevent since it should always be accessed via the Python interface provided by :mod:`a5py`.

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

.. _Scripts:

Terminal scripts
================

Some useful scripts found in `./bin` folder are accessible from terminal if ``a5py`` was installed succesfully.
These scripts modify the contents of the HDF5 file and they are listed here.

.. list-table:: Command line scripts
   :widths: 25 75
   :header-rows: 0

   * - ``a5editoptions``
     - Edit currently active simulation options (or create default options if none exists).
       Options are edited with a text editor and set as active after saving.
   * - ``a5manage``
     - Manage data groups by setting them as active, removing them, or copying them to another file.
       Inputs that have been used by a run cannot be removed unless that run is removed first.
   * - ``a5gui``
     - Open GUI where you can manage your data and view results and inputs.
   * - ``a5combine``
     - Combine simulation results if: a run was split into several processes without MPI (this creates multiple output files), several simulations were run in parallel and stored in different files, an existing simulation was continued.
   * - ``a5upgrade``
     - Makes a copy "<originalname>_<currentversion>.h5" of the file if it was made with an older version of ASCOT5.
       The aim is to maintain backwards compatibility with this script.
       However, it is not always possible to convert the results into a new format and in those cases rerunning the simulation with updated inputs is unfortunately necessary.

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
