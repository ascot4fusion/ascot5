import numpy as np
import unyt

import pytest

from a5py import Ascot, SimulationOptions
from a5py import plotting as a5plt


def test_no_drifts(ascot, inspect, plot):
    if not inspect:
        parameters = SimulationOptions(
                simulation_mode=1,
                enable_adaptive=False,
                timestep=1e-11,
                enable_orbit_following=True,
                collect_orbit=True,
                number_of_points_per_marker=202,
                interval=1e-11,
                activate_simulation_time_limits=True,
                max_mileage=2e-9,
            )

        ascot.data.create_guidingcentermarker(
            species="electron",
            ids=np.array([1, 2]),
            r=5*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=100e6*unyt.eV,
            pitch=0.5,
            gyroangle=0.0*unyt.rad,
            charge=np.array([1, -1])*unyt.e,
            )

        ascot.data.create_bfieldcartesian(
            bxyz=np.array([5., 0., 0.])*unyt.T,
            jacobian=np.full((3, 3), 0.)*unyt.T/unyt.m,
            axisrz=np.array([1., 0.])*unyt.m,
            rhoval=1.,
            )

        ascot.data.create_efieldcartesian(
            exyz=np.array([0, 0, 0])*unyt.V/unyt.m,
            )

        ascot.data.bfield.active.stage()
        ascot.data.efield.active.stage()
        ascot.data.marker.active.stage()
        run = ascot.simulate(params=parameters)

    x1, y1, z1 = run.getorbit("x", "y", "z", mode="particle", filter=[1])
    x2, y2, z2 = run.getorbit("x", "y", "z", mode="particle", filter=[2])
    if plot:
        a5plt.plot(y1,z1)
        a5plt.plot(y2,z2)
        a5plt.show()


@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_exb_drift(method):
    a5 = Ascot()

    parameters = SimulationOptions(
        simulation_mode=1,
        enable_adaptive=False,
        enable_orbit_following=True,
        collect_orbit=True,
        activate_simulation_time_limits=True,
        max_mileage=1e-7,
        )

    if method == "go":
        parameters.simulation.simulation_mode = 1
        parameters.simulation.timestep=1e-10
        parameters.orbit.number_of_points_per_marker=10002
        parameters.orbit.interval=1e-10
    else:
        parameters.simulation.simulation_mode = 2
        parameters.orbit.number_of_points_per_marker=1002
        parameters.orbit.interval=1e-9
        parameters.simulation.timestep = 1e-9
        if method == "gc-fixedstep":
            parameters.simulation.enable_adaptive = False

    a5.data.create_guidingcentermarker(
        species="electron",
        ids=np.array([1, 2]),
        r=5*unyt.m,
        z=0.0*unyt.m,
        phi=0.0*unyt.deg,
        ekin=100e6*unyt.eV,
        pitch=0.5,
        gyroangle=0.0*unyt.rad,
        charge=np.array([1, -1])*unyt.e,
        )

    a5.data.create_bfieldcartesian(
        bxyz=np.array([5., 0., 0.])*unyt.T,
        jacobian=np.full((3, 3), 0.)*unyt.T/unyt.m,
        axisrz=np.array([1., 0.])*unyt.m,
        rhoval=1.,
        )

    a5.data.create_efieldcartesian(
        exyz=np.array([0, 1e6, 0])*unyt.V/unyt.m,
        )

    run = a5.simulate(params=parameters)

    x1, y1, z1 = run.getorbit("x", "y", "z", mode="particle", filter=[1])
    x2, y2, z2 = run.getorbit("x", "y", "z", mode="particle", filter=[2])
    import matplotlib.pyplot as plt
    plt.plot(y1,z1)
    plt.plot(y2,z2)
    plt.show()

def test_gradb_drift():
    a5 = Ascot()
    parameters = SimulationOptions(
        simulation_mode=1,
        enable_adaptive=False,
        use_explicit_fixedstep=True,
        explicit_fixedstep=1e-10,
        enable_orbit_following=True,
        collect_orbit=True,
        number_of_points_per_marker=10002,
        interval=1e-11,
        activate_simulation_time_limits=True,
        max_mileage=1e-7,
        activate_real_time_limit=True,
        max_real_time=5.0,
    )

    a5.data.create_guidingcentermarker(
        species="electron",
        ids=np.array([1, 2]),
        r=5*unyt.m,
        z=0.0*unyt.m,
        phi=0.0*unyt.deg,
        ekin=100e6*unyt.eV,
        pitch=0.5,
        gyroangle=0.0*unyt.rad,
        charge=np.array([1, -1])*unyt.e,
        )

    jacobian = np.array(
        [
            [0., 0., 0.1],
            [0., 0., 0.],
            [0., 0., 0.],
        ]
        )*unyt.T/unyt.m
    a5.data.create_bfieldcartesian(
        bxyz=np.array([5., 0., 0.])*unyt.T,
        jacobian=jacobian,
        axisrz=np.array([1., 0.])*unyt.m,
        rhoval=1.,
        )

    a5.data.create_efieldcartesian(
        exyz=np.array([0, 0, 0])*unyt.V/unyt.m,
        )

    run = a5.simulate(params=parameters)

    #print(run.getorbit("x", "y", "z", mode="guidingcenter", filter=[1]))
    x1, y1, z1 = run.getorbit("x", "y", "z", mode="particle", filter=[1])
    x2, y2, z2 = run.getorbit("x", "y", "z", mode="particle", filter=[2])
    #print(z1.size, z2.size)
    import matplotlib.pyplot as plt
    plt.plot(y1,z1)
    plt.plot(y2,z2)
    plt.show()
