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

::

    data
    ├── bfield               # Magnetic field inputs
    │   ├── B_2DS_7027705680 # Some 2D magnetic field data
    │   ├── B_3DS_0890178582 # Some other 2D magnetic field data
    │   └── ...              # Some other possible magnetic field data
    │
    ├── efield               # Electric field inputs
    │   └── ...
    │
    ├── ...                  # Other inputs (wall, plasma, etc.)
    │
    ├── run_892758002        # Results of a single simulation
    ├── run_992765110        # Results of another simulation
    └── ...                  # Other possible results

The input data is divided into separate *parent groups*: one for the magnetic field inputs, one for the electric field inputs and so on.
Each parent group can have multiple children which contain the actual input data.
For example, one may have several plasma inputs that vary in density for a parameter scan or one can have a separate 2D and 3D magnetic fields for comparison.
One of the children is always marked as **active** meaning that input will be used in the next simulation.
The active input is used by default when interpolating the input data via the Python interface in pre- and postprocessing.

As for the results, one group is always marked as active (if any results exist) which by default is the result of the most recent simulation.
Each input and result has a ten number string that is randomly generated when that data is written to the file.
This **QID** (quasi unique identifier) is used to separate different inputs and results from one another.
Each input and result also contains the date at which they were created, what type they are (e.g. 2D vs 3D tokamak magnetic field), and an optional description given by the user.
It is strongly recommended to document your work using the description field which can also be used to conveniently access the data in postprocessing.

Data access
===========

.. admonition:: Summary

   - Data is accessed, generated, and processed via :class:`.Ascot` instance.
   - The attribute ``Ascot.data``, which is an instance of :class:`.Ascot5IO`, is used to traverse the contents of the file.
   - Inputs and results can be referred either by their full name, QID, or user-defined tag that is the first word on the description on all caps.
   - Inputs and results are instances of :class:`.DataGroup` and contain methods to set the description or activate or remove the group from the file.
   - All changes, e.g. setting group as active, persist between sessions as the information is stored on the file.

Instances of :class:`.Ascot` class are used to open ``ascot.h5`` or create an empty file.

.. code-block:: python

   from a5py import Ascot
   a5 = Ascot("ascot.h5", create=True) # Creates a new file
   a5 = Ascot("ascot.h5")              # Opens an existing file

No actual data is read during the initialization, so these objects are very light-weight and new instances can be created at will
The contents of the file are accessed via ``Ascot.data``, which is an instance of :class:`.Ascot5IO`, and it provides a "treeview" of the file which was illustrated in the previous section.

.. warning::
   :class:`.Ascot` does not monitor if the data has changed on the disk.
   If the data was changed by another instance or by running a simulation, the object must be reinitialized.
   No reinitialization is needed if the object itself changes the contents of the file as those are monitored.

Inputs are divided to groups depending on what type of input data they provide.
All magnetic field inputs are located in the ``bfield`` group, all electric field inputs in the ``efield`` group, and so on.
All different *input parent* groups are instances of :class:`.InputNode` class.

One of the inputs in each group is always flagged active in the file, meaning it will be used as an input the next time a simulation is run.
For example, the currently active magnetic field input would be

.. code-block:: python

   a5.data.bfield.active # <- Currently active magnetic field input

The magnetic field object here, like the results and input data in general, are stored as objects that inherit from the :class:`.DataGroup` class.

Similarly, results groups are instances of :class:`.ResultNode` class and one of them is always flagged active.
Whenever a simulation is ran, the newly created group is always set as active.
Simulation results are located at the top level of the hierarchy, e.g. the :class:`.DataGroup` corresponding to the active run is

.. code-block:: python

   a5.data.active # <- Currently active run (most recent by default)

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
The first word (omitting special characters) of the description is used to form a **tag** for that group which can be used to make accessing the data more convenient:

.. code-block:: python

   a5.data.get_date() # Returns date when this group (run) was created.
   a5.data.set_desc("Alpha particle simulation")
   a5.data.ALPHA      # Refers to the group on the previous line.

If multiple groups have same tag, they are separated by a running index which is added to the tag.
The group can be set active or removed from the file with the corresponding methods:

.. code-block:: python

   a5.data.bfield.MYFIELD.activate() # Sets group with tag "MYFIELD" as active
   a5.data.bfield.MYFIELD.destroy()  # Removes the group from the file.
   a5.data.bfield.destroy()          # Removes all magnetic field inputs
   a5.data.destroy()                 # Removes all results

Finally, the treeview can be accessed both attribute-like and dictionary-like.

.. code-block:: python

   # These are all equivalent references to a group
   a5.data.bfield.q0123456789
   a5.data["bfield"]["q0123456789"]
   a5.data.bfield["TAG"]

.. rubric:: References

.. autosummary::
   :nosignatures:

   ~a5py.ascot5io.Ascot5IO
   ~a5py.ascot5io.RunGroup
   ~a5py.ascot5io.InputGroup
   ~a5py.ascot5io.coreio.treedata.DataGroup

.. _Scripts:

Terminal scripts
================

Some useful scripts can be found in `./bin` folder or accessed directly from terminal if ``a5py`` was installed succesfully.
These scripts modify or disbplay the contents of the HDF5 file and they are listed here.

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
