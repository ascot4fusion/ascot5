"""Test magnetic field module."""
import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.physlib.formulas import cart2pol, cart2pol_vec
from a5py.engine.interpolate import evaluate
from a5py.templates import (
    PremadeMagneticField, TransmuteAnalyticalBfieldToSplines,
)
from a5py.physlib import aeq

from .conf import test_data, create_input

@pytest.fixture
def ascot():
    return Ascot()


def test_bfieldcartesian_init(ascot):
    obj = ascot.data.create_bfieldcartesian(
        bxyz=test_data["BfieldCartesian"]["bxyz"],
        jacobian=test_data["BfieldCartesian"]["jacobian"],
        axisrz=test_data["BfieldCartesian"]["axisrz"],
        rhoval=test_data["BfieldCartesian"]["rhoval"],
        )
    assert obj.psival.v == obj.rhoval.v


def test_interpolate_cartesian():
    jac = np.array([
        [1., 2., 3.],
        [4., 5., 6.],
        [7., 8., 9.],
        ]) * unyt.T/unyt.m
    obj = create_input("BfieldCartesian", jacobian=jac)
    x, y, z = 2.*unyt.m, 0.*unyt.m, 0.*unyt.m
    r, phi, z = cart2pol(x, y, z)
    br0, bphi0, bz0, rho, psi = evaluate(
        r, phi, z, 0.*unyt.s, "br", "bphi", "bz", "rho", "psi", bfield=obj
        )

    jac = obj.jacobian
    bx = obj.bxyz[0] + jac[0,0]*x + jac[0,1]*y + jac[0,2]*z
    by = obj.bxyz[1] + jac[1,0]*x + jac[1,1]*y + jac[1,2]*z
    bz = obj.bxyz[2] + jac[2,0]*x + jac[2,1]*y + jac[2,2]*z

    br, bphi, bz = cart2pol_vec(x, bx, y, by, z, bz)
    assert np.isclose(br, br0)
    assert np.isclose(bz, bz0)
    assert np.isclose(bphi, bphi0)
    assert np.isclose(rho, obj.rhoval)
    assert np.isclose(psi, obj.psival)

    jac0 = evaluate(
        r, phi, z, 0.*unyt.s, "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz", bfield=obj
        )
    dr = 1.e-1 * unyt.m
    br, bphi, bz = evaluate(
        r+dr, phi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[0] * dr)
    assert np.isclose(bphi - bphi0, jac0[3] * dr)
    assert np.isclose(bz - bz0, jac0[6] * dr)

    dz = 1.e-1 * unyt.m
    br, bphi, bz = evaluate(
        r, phi, z+dz, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[2] * dz)
    assert np.isclose(bphi - bphi0, jac0[5] * dz)
    assert np.isclose(bz - bz0, jac0[8] * dz)

    dphi = 1.e-3 * unyt.deg
    br, bphi, bz = evaluate(
        r, phi+dphi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    dphi = dphi.to("rad").v
    assert np.isclose(br - br0, jac0[1] * dphi)
    assert np.isclose(bphi - bphi0, jac0[4] * dphi)
    assert np.isclose(bz - bz0, jac0[7] * dphi)


def test_interpolate_analytical(ascot):
    template = PremadeMagneticField(ascot, field="iter-baseline")
    obj = template.create_input()

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n))*unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n))*unyt.m

    rho0, psi0 = evaluate(
        r, 0*unyt.deg, z, 0.*unyt.s, "rho", "psi", bfield=obj
        )
    r0 = obj.rmajor
    psi = aeq.psi(r/r0, z/r0, obj.coefficients[:-1], obj.coefficients[-1]) * obj.psiscaling
    assert np.isclose(psi, psi0).all()
    rho = np.sqrt( (psi - obj.psilimits[0]) / (obj.psilimits[1] - obj.psilimits[0]) )
    assert np.isclose(rho, rho0).all()

    r0, z0, psi0 = r[0], z[0], psi[0]
    dr, dz = 1.e-4 * unyt.m, 1.e-4 * unyt.m
    psi = evaluate(
        unyt.unyt_array([r0+dr, r0]), 0*unyt.deg, unyt.unyt_array([z0, z0+dz]),
        0.*unyt.s, "psi", bfield=obj
        )

    bz = (unyt.rad*(psi[0] - psi0) / (r0*dr)).to("T")
    br = -(unyt.rad*(psi[1] - psi0) / (r0*dz)).to("T")
    br0, bz0 = evaluate(r0, 0*unyt.deg, z0, 0.*unyt.s, "br", "bz", bfield=obj)
    assert np.isclose(br, br0, atol=1e-3)
    assert np.isclose(bz, bz0, atol=1e-3)

    r, z, phi = r0, z0, 0.0*unyt.deg
    br0, bphi0, bz0 = evaluate(
        r, phi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )

    jac0 = evaluate(
        r, phi, z, 0.*unyt.s, "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz", bfield=obj
        )
    dr = 1.e-4 * unyt.m
    br, bphi, bz = evaluate(
        r+dr, phi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[0] * dr)
    assert np.isclose(bphi - bphi0, jac0[3] * dr)
    assert np.isclose(bz - bz0, jac0[6] * dr)

    dz = 1.e-4 * unyt.m
    br, bphi, bz = evaluate(
        r, phi, z+dz, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[2] * dz)
    assert np.isclose(bphi - bphi0, jac0[5] * dz)
    assert np.isclose(bz - bz0, jac0[8] * dz)

    dphi = 1.e-4 * unyt.deg
    br, bphi, bz = evaluate(
        r, phi+dphi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    dphi = dphi.to("rad").v
    assert np.isclose(br - br0, jac0[1] * dphi)
    assert np.isclose(bphi - bphi0, jac0[4] * dphi)
    assert np.isclose(bz - bz0, jac0[7] * dphi)


def test_interpolate_spline2d(ascot):
    template = PremadeMagneticField(ascot, field="iter-baseline")
    correct = template.create_input()

    template = TransmuteAnalyticalBfieldToSplines(
        ascot, field=correct, nr=400, rlim=unyt.unyt_array([4, 8.2], "m"),
        nz=800, zlim=unyt.unyt_array([-4, 4], "m"))
    testdata = template.create_input()

    n = 100
    #r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    #z = (-4 + (4 - -4) * np.random.rand(n)) * unyt.m
    r, z = 6.5 * unyt.m, 0.0*unyt.m

    br0, bphi0, bz0, rho0, psi0 = evaluate(
        r, 0*unyt.deg, z, 0.*unyt.s, "br", "bphi", "bz", "rho", "psi", bfield=correct
        )

    br, bphi, bz, rho, psi = evaluate(
        r, 0*unyt.deg, z, 0.*unyt.s, "br", "bphi", "bz", "rho", "psi", bfield=testdata
        )

    assert np.isclose(br, br0, rtol=1e-1).all()
    assert np.isclose(bz, bz0, rtol=1e-1).all()
    assert np.isclose(bphi, bphi0, rtol=1e-1).all()
    assert np.isclose(rho, rho0, rtol=1e-1).all()
    assert np.isclose(psi, psi0, rtol=1e-1).all()

    jac0 = evaluate(
        r, 0*unyt.deg, z, 0.*unyt.s, "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz", bfield=correct
        )

    jac = evaluate(
        r, 0*unyt.deg, z, 0.*unyt.s, "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz", bfield=testdata
        )

    for qnt in zip(jac, jac0):
        assert np.isclose(qnt[0], qnt[1], rtol=1e-1).all()
