"""Verify that the code reproduces gyromotion and guiding center drifts
accurately.

The tests in this file are done with energetic electrons to verify that the
orbit-following physics work in the relativistic regime.
"""
import numpy as np
import unyt

import pytest

from a5py import SimulationOptions, physlib
from a5py import plotting as a5plt

from helpers import assert_isclose


def test_drifts_nodrifts(ascot, inspect, plot):
    """Verify that the code reproduces correct gyroradius and gyrofrequency in
    a uniform magnetic field.
    """
    if not inspect:
        parameters = SimulationOptions.from_dict(
            simulation={
                "mode": "gyro-orbit",
                "timestep": 1e-11,
            },
            physics={
                "enable_orbit_following": True,
            },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 2e-9,
            },
            orbit={
                "collect": "interval",
                "buffer_size": 202,
                "interval": 1e-11,
            },
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

        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    m, q, pitch, ekin, x0, y0 = run.getstate("mass", "charge", "pitch", "ekin", "y", "z")
    x, y, time = run.getorbit("y", "z", "time", filter=[1])

    gyro_radius, gyro_frequency = {}, {}
    bnorm = np.sqrt(np.sum(run.bfield.bxyz**2))
    gamma = physlib.Quantity.registry["gamma"].compute(mass=m, ekin=ekin)
    gyro_radius["analytical"] = (
        physlib.Quantity.registry["gyroradius"].compute(
            mass=m, charge=q, gamma=gamma, pitch=pitch, bnorm=bnorm,
            ).to("m")
    )
    gyro_frequency["analytical"] = (
        physlib.Quantity.registry["gyrofrequency"].compute(
            mass=m, charge=q, bnorm=bnorm, gamma=gamma,
            ).to("rad/s")
    )
    gyro_radius["numerical"] = (
        np.mean( np.sqrt( (x - x0[0])**2 + (y - y0[0])**2 ) )
    )
    gyro_frequency["numerical"] = (
        np.sum( np.sqrt( np.diff(x)**2 + np.diff(y)**2 ) )
        * unyt.rad / (gyro_radius["numerical"] * time[-1])
    )

    if plot:
        _, ax = a5plt.subplots()

        ax.set_ylabel("z [m]")
        ax.set_xlabel("y [m]")

        ax.set_aspect("equal", adjustable="box")
        orbx = x - x0[0]
        orby = y - y0[0]

        anax = gyro_radius["numerical"] * np.sin(np.linspace(0, 2*np.pi, 360))
        anay = gyro_radius["numerical"] * np.cos(np.linspace(0, 2*np.pi, 360))
        fangle = (
            np.arctan2(orby[0], orbx[0]) * unyt.rad + gyro_frequency["analytical"][0] * time[-1]
        )
        fx = gyro_radius["numerical"] * np.cos(-fangle) # Minus sign because this we are
        fy = gyro_radius["numerical"] * np.sin(-fangle) # in LHS coordinate system.

        ax.plot(orbx, orby, linewidth=3)
        ax.plot(anax, anay, linestyle="--", color="black", alpha=0.7)
        ax.scatter(orbx[-1], orby[-1], marker="o", color="C0", zorder=3)
        ax.scatter(fx, fy, marker="x", color="black", zorder=4)
        a5plt.show()

    assert_isclose("Gyroradius", gyro_radius["numerical"], gyro_radius["analytical"][0], 1e-3)
    assert_isclose("Gyrofrequency", gyro_frequency["numerical"], gyro_frequency["analytical"][0], 1e7)


@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_drifts_exb(ascot, method, inspect, plot):

    if not inspect:
        parameters = SimulationOptions.from_dict(
            physics={
                "enable_orbit_following": True,
            },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 1e-7,
            },
            orbit={
                "collect": "interval",
            },
        )

        if method == "go":
            parameters.simulation.mode = "gyro-orbit"
            parameters.simulation.timestep = 1e-10
            parameters.orbit.buffer_size = 10002
            parameters.orbit.interval = 1e-10
        else:
            parameters.simulation.mode = "guiding-center"
            parameters.simulation.timestep = 1e-9
            parameters.orbit.buffer_size = 1002
            parameters.orbit.interval = 1e-9
            if method == "gc-fixedstep":
                parameters.simulation.enable_adaptive = False

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
            exyz=np.array([0, 1e6, 0])*unyt.V/unyt.m,
            )

        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    yi = run.getstate("z", state="ini")
    yf, deltat = run.getstate("z", "time", state="end")
    x1, y1 = run.getorbit("y", "z", filter=[1])
    x2, y2 = run.getorbit("y", "z", filter=[2])

    # Analytical values
    bvec = run.bfield.bxyz.ravel()
    evec = run.efield.exyz.ravel()
    v_ExB = np.cross(evec, bvec) / np.inner(bvec, bvec)
    time = run.options.endconditions.max_mileage * unyt.s

    # Numerical values
    vnum_ExB = (yf - yi) / deltat

    if plot:

        def plot_velocity(x, y, dt, text):
            ini = np.array([x.v + .01, y.v])
            end = np.array([x.v + .01, (y + v_ExB[2]*dt).v])
            ax.annotate("", xytext=ini, xy=end,
                        arrowprops={"arrowstyle":"->", "color":"black"})
            ax.plot([ini[0], end[0]], [ini[1], end[1]], color="black")
            ax.annotate(text, ini + .005)

        _, ax = a5plt.subplots()
        ax.plot(x1, y1)
        ax.plot(x2, y2)

        plot_velocity(x1[0], y1[0], deltat[0], r"$e^-$")
        plot_velocity(x2[0], y2[0], deltat[1], r"$e^+$")

        a5plt.show()

    assert_isclose("E x B drift (positron)", vnum_ExB[2], v_ExB[0], 1e-3)
    assert_isclose("E x B drift (electron)", vnum_ExB[2], v_ExB[1], 1e-3)


@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_drifts_gradb(ascot, method, inspect, plot):
    if not inspect:
        parameters = SimulationOptions.from_dict(
            physics={
                "enable_orbit_following": True,
            },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 1e-7,
            },
            orbit={
                "collect": "interval",
            },
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

        jacobian = np.array(
            [
                [0., 0., 0.1],
                [0., 0., 0.],
                [0., 0., 0.],
            ]
            )*unyt.T/unyt.m
        ascot.data.create_bfieldcartesian(
            bxyz=np.array([5., 0., 0.])*unyt.T,
            jacobian=jacobian,
            axisrz=np.array([1., 0.])*unyt.m,
            rhoval=1.,
            )

        ascot.data.create_efieldcartesian(
            exyz=np.array([0, 0, 0])*unyt.V/unyt.m,
            )

        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    pitch, ekin, q, m = run.getstate("pitch", "ekin", "charge", "mass", filter=[1])
    xi = run.getstate("y", state="ini")
    xf, deltat = run.getstate("y", "time", state="end")
    x1, y1 = run.getorbit("y", "z", filter=[1])
    x2, y2 = run.getorbit("y", "z", filter=[2])

    # # Analytical values
    bvec = run.bfield.bxyz.ravel()
    gradb = run.bfield.jacobian[0,:]
    bnorm = np.sqrt(np.sum(bvec**2))
    gamma = physlib.Quantity.registry["gamma"].compute(mass=m, ekin=ekin)
    vnorm = physlib.Quantity.registry["vnorm"].compute(gamma=gamma)
    v_gradb = (
        gamma * m * ( (1.0 - pitch**2) * vnorm**2 / ( q * bnorm ) )
            * 0.5 * np.cross(bvec, gradb) / bnorm**2
    ).to("m/s")
    time = run.options.endconditions.max_mileage * unyt.s

    # Numerical values
    vnum_gradb = (xf - xi) / deltat

    if plot:
        def plot_velocity(x, y, dt, text):
            ini = np.array([x.v + .01, y.v])
            end = np.array([x.v + .01, (y + v_gradb[1]*dt).v])
            ax.annotate("", xytext=ini, xy=end,
                        arrowprops={"arrowstyle":"->", "color":"black"})
            ax.plot([ini[0], end[0]], [ini[1], end[1]], color="black")
            ax.annotate(text, ini + .005)

        _, ax = a5plt.subplots()
        ax.plot(x1, y1)
        ax.plot(x2, y2)

        plot_velocity(x1[0], y1[0], deltat[0], r"$e^-$")
        plot_velocity(x2[0], y2[0], deltat[1], r"$e^+$")
        a5plt.show()

    assert_isclose("E x B drift (positron)", vnum_gradb[0], v_gradb[1], 1e-3)
    assert_isclose("E x B drift (electron)", vnum_gradb[1], v_gradb[1], 1e-3)
