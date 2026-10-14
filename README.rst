ASCOT5
======

`ASCOT5 <https://arxiv.org/abs/1908.02482/>`_ (or simply *ASCOT*) is a test-particle orbit-following software for computing particle orbits in three-dimensional magnetic-field geometries.
Its output includes particle orbits, phase-space distributions, and wall loads.
ASCOT is frequently used to study fast ions, trace impurities, and runaway electrons in tokamaks, stellarators, and other magnetically confined plasmas.
Particle orbits can be followed either fully, including the gyro-motion, or using a reduced model in which only the guiding-center trajectory is followed.

The user interface and API are implemented in Python, while the computationally intensive parts of the code are written in C.
ASCOT is extensively parallelized and optimized for both CPUs and GPUs.
Simulation data can be stored in HDF5 format.
ASCOT can also be directly integrated with other modelling frameworks, such as `IMAS <https://doi.org/10.1088/0029-5515/55/12/123006>`_.

ASCOT is open source and `the repository <https://github.com/ascot4fusion/ascot5>`_ is maintained by `VTT Technical Research Centre of Finland <https://www.vttresearch.com/en/ourservices/future-nuclear#fusion-reactor-research>`_.
For contact: Konsta Särkimäki and Antti Snicker (firstname.lastname@vtt.fi).

The team at VTT also offers technical support and modeling work using ASCOT.

Using ASCOT
===========

Installing and using ASCOT is extensively covered in the `documentation <https://ascot4fusion.github.io/ascot5/>`_.
For a quick start, install the code using conda

.. code-block:: bash

    git clone https://github.com/ascot4fusion/ascot5.git
    cd ascot5
    conda env create -f environment.yaml
    conda activate ascot-env
    make ascot
    pip install .

and then you can try the code by simulating a particle in a uniform magnetic field:

.. code-block:: python

    import unyt
    import numpy as np
    from a5py import Ascot, SimulationOptions

    ascot = Ascot()
    ascot.data.create_bfieldcartesian(
        bxyz=np.array([0., 1., 0])*unyt.T,
        jacobian=np.zeros((3,3))*unyt.T/unyt.m,
        axisrz=[1., 0.]*unyt.m,
        rhoval=0,
        )
    ascot.data.create_efieldcartesian(
        exyz=np.zeros((3,))*unyt.V/unyt.m,
        )
    ascot.data.create_guidingcentermarker(
        r=1.*unyt.m,
        z=0.*unyt.m,
        phi=0.*unyt.deg,
        ekin=1e6*unyt.eV,
        pitch=0.5,
        gyroangle=0.*unyt.rad,
        species="H",
        )

    parameters = SimulationOptions.from_dict(
        simulation={
            "mode": "gyro-orbit",
            "timestep": 1e-10,
        },
        physics={
            "enable_orbit_following": True,
        },
        endconditions={
            "activate_simulation_time_limits": True,
            "max_mileage": 1e-6,
        },
        orbit={
            "collect": "interval",
            "buffer_size": int(1e-6/1e-10),
            "interval": 1e-10,
        },
    )

    run = ascot.simulate(params=parameters)
    run.plotorbit_trajectory("y", "x")

Licence
=======

The ASCOT5 and associated programs are distributed under the terms of the GNU Lesser General Public License (LGPL).
Please see the files COPYING and COPYING.LESSER for more information.
