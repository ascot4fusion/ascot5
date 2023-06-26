==========
Python API
==========

``a5py``
========

.. automodule:: a5py

.. toctree::
   :hidden:

   a5py.ascot5io
   a5py.ascotpy
   a5py.templates
   a5py.physlib
   a5py.routines
   a5py.gui

.. autosummary::

   a5py.ascot5io
   a5py.ascotpy
   a5py.templates
   a5py.physlib
   a5py.routines
   a5py.gui

``exceptions``
===================

.. automodule:: a5py.exceptions
   :members:

``Ascot``
==============

.. autoclass:: a5py.Ascot

Input initialization
********************

.. autosummary::

   a5py.Ascot.input_init
   a5py.Ascot.input_free

Live simulations
****************

.. autosummary::

   a5py.Ascot.input_init
   a5py.Ascot.input_free

API
***

.. automethod:: a5py.Ascot.input_init
.. automethod:: a5py.Ascot.input_free
.. automethod:: a5py.Ascot.simulation_initinputs
.. automethod:: a5py.Ascot.simulation_run
