.. Main page of the user documentation

User guide
==========

.. image:: ../_static/promo.svg

----------------------------------

`ASCOT5 <https://arxiv.org/abs/1908.02482/>`_ (or simply *ASCOT*) is a test-particle orbit-following software for computing particle orbits in three-dimensional magnetic-field geometries.
Its output includes particle orbits, phase-space distributions, and wall loads.
ASCOT is frequently used to study fast ions, trace impurities, and runaway electrons in tokamaks, stellarators, and other magnetically confined plasmas.
Particle orbits can be followed either fully, including the gyro-motion, or using a reduced model in which only the guiding-center trajectory is followed.

The user interface and API are implemented in Python, while the computationally intensive parts of the code are written in C.
ASCOT is extensively parallelized and optimized for both CPUs and GPUs.
Simulation data can be stored in HDF5 format.
ASCOT can also be directly integrated with other modelling frameworks, such as `IMAS <https://doi.org/10.1088/0029-5515/55/12/123006>`_.

ASCOT is open source and `the repository <https://github.com/ascot4fusion/ascot5>`_ is maintained by `VTT Technical Research Centre of Finland <https://www.vttresearch.com/en/ourservices/fusion-energy-research-and-services>`_.
The team at VTT also offers technical support and modeling work using ASCOT.

For contact (firstname.lastname@vtt.fi):

.. grid:: 1 2 3 4

    .. grid-item::

        .. image:: ../_static/profile_konsta.png
            :height: 200px
            :align: center

        **Konsta Särkimäki**

        Lead Developer

    .. grid-item::

        .. image:: ../_static/profile_antti.png
            :height: 200px
            :align: center

        **Antti Snicker**

        Group Leader

.. note::

    If you have done simulations where you compare ASCOT to experimental results or output of other similar software, please contact the developers!
    We would like to maintain a record of those.

Getting help
************

For questions related to the code or its physics, please join our `Slack <http://ascot-workspace.slack.com/>`_.
We also hold ASCOT meetings every two weeks.
These meetings are primarily intended for discussing ongoing or planned studies involving ASCOT, and everyone is welcome to join and discuss their work with other ASCOT users and experts.
In addition, office hours ("ASCOT Garage") are held on most Wednesdays for questions related to using the code.
The schedules and connection information for both the ASCOT meetings and the Garage are available in Slack, in the pinned post in the ascot5 channel.

License
*******

This software is distributed under the `LGPLv3 <https://www.gnu.org/licenses/lgpl-3.0.html>`_ license.

.. toctree::
    :maxdepth: 2
    :caption: Getting started

    setup

.. toctree::
    :maxdepth: 2
    :caption: Using ASCOT5

    tutorials
    physics
    testing
    publications

.. toctree::
    :maxdepth: 2
    :caption: API reference:

    dataform
    inputgen
    postproc
