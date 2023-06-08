===========================
Data access and preparation
===========================

Simulation inputs and options, as well as the outputs once the simulation is complete, are stored in a single HDF5 file.

The ASCOT5 HDF5 file, `ascot.h5` from here on out, contains all data required to run a simulation, and so the results and the inputs that produced them are always kept together.
The file is designed to hold multiple inputs (even of same type) and results of multiple simulations to encourage one to use a single file for the whole project or study.

The exact structure of the HDF5 file is irrelevant since it should always be accessed via the Python interface provided by `a5py`.

Data access
===========

Instances of :class:`.Ascot` class are used to open ``ascot.h5``.
The contents of the file are accessed via `Ascot.data`, which is an instance of :class:`.Ascot5IO`, and it provides a "treeview" of the file.
For example, the currently *active* magnetic field input is

The magnetic field object here, like the results and input data in general, are stored as objects that inherit from the :class:`.DataGroup` class.
These objects don't contain any data, but provide an interface to access the data on the disk.
This means that initializing the `Ascot` object is a very light-weight process since no actual data is read during the initialization.

Simulation results are located at the top level of the hierarchy, e.g. the :class:`.DataGroup` corresponding to the *active* run is

Inputs are divided on groups depending on what type of input data they provide. All magnetic field inputs are located under the `bfield` group, all electric field inputs under the `efield`, and so on.
All different *input parent* groups are listed in `INPUGROUPS` and they are instances of :class:`.InputNode` class.
One of the inputs in each group is always flagged *active*, meaning it will be used as an input the next time a simulation is run.
Similarly one of the simulation results is always flagged active (by default it is the latest) and they are instances of :class:`.ResultNode` class.
The active inputs are also used in post-processing and in GUI, so pay attention to what groups are active at the moment.
The group can be set active with

The active flag is stored in the HDF5 file so it persists between Python sessions.

Each group has a QID which is a 10-number identifier to separate that group from others and it is automatically created.
For humans there is a description attribute that one can use to document contents of `ascot.h5`.

The meta data is accessed and modified using the methods found in :class:`.DataGroup` where one can also find methods to remove data groups.
See the tutorial and API below for detailed description of data access.

.. autoclass:: a5py.ascot5io.Ascot5IO
  :members:
  :inherited-members:
  
.. autoclass:: a5py.ascot5io._iohelpers.treeview.InputNode
  :members:
  :inherited-members:
  
.. autoclass:: a5py.ascot5io._iohelpers.treeview.ResultNode
  :members:
  :inherited-members:
  
.. autoclass:: a5py.ascot5io._iohelpers.treedata.DataGroup
  :members:
  :inherited-members:


Input generation
================




Post-processing
===============


.. autoclass:: a5py.Ascot
  :members:
  :inherited-members:

Interactive simulations
=======================


Exceptions
==========
