"""Test magnetic field module."""
import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.physlib.formulas import cart2pol, cart2pol_vec
from a5py.engine.interpolate import evaluate

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


def test_interpolate(ascot):
    jac = np.array([
        [1., 2., 3.],
        [4., 5., 6.],
        [7., 8., 9.],
        ]) * unyt.T/unyt.m
    obj = create_input("BfieldCartesian", jacobian=jac)
    x, y, z = 1.*unyt.m, 2.*unyt.m, 3.*unyt.m
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
