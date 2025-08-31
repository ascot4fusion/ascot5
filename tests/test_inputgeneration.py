"""Test input generation."""
import os

import pytest
import numpy as np

from a5py.data import AscotData

from unyt import dimensionless as nodim, T, m, Wb, deg, s, V


variants = [
        ("bfieldcartesian", 1),
        ("bfieldanalytical", 5),
        ("bfield2d", 2),
        ("bfield3d", 2),
        ("bfieldstellarator", 3),
        ("efieldcartesian", 0),
        ("efieldradialpotential", 0),
        ("plasma1d", 4),
        ("plasma1ddynamic", 4),
        ("neutral1d", 0),
        ("neutral3d", 0),
        ("mhdstationary", 3),
        ("mhddynamic", 3),
        #("boozermap", 0),
        #("atomicdata", 0),
        ("wall2d", 2),
        ("wall3d", 2),
        #("nbi", 0),
        #("particlemarker", 0),
        #("guidingcentermarker", 0),
        #("fieldlinemarker", 0),
    ]
"""List of tuples of variants and number of optional parameters in each."""

pytestmark = pytest.mark.parametrize(
    "variant, nopt", variants,
    ids=[name[0] for name in variants],
)


@pytest.fixture(name="datamanager")
def fixture_datamanager():
    yield AscotData(hdf5file=("testascot.h5", False))
    os.unlink("testascot.h5")


@pytest.mark.parametrize(
    "store_hdf5",
    [True, False],
    ids=["store_hdf5", "dont_store_hdf5"],
    )
@pytest.mark.parametrize(
    "include_opt",
    [True, False],
    ids=["include_optional", "exclude_optional"],
    )
def test_create(
    datamanager, variant, store_hdf5, nopt, include_opt
    ):
    """Test input creation."""
    create_input = getattr(datamanager, "create_" + variant)
    obj = create_input(dryrun=True)
    properties = obj.export()

    while nopt and not include_opt:
        nopt -= 1
        properties.popitem()
    obj = create_input(**properties, store_hdf5=store_hdf5)
    for property, value in properties.items():
        if hasattr(value, "units"):
            assert np.allclose(
                getattr(obj, property), value, atol=0*value.units
                )
        elif isinstance(value, dict):
            assert getattr(obj, property) == value
        elif isinstance(value[0], str):
            assert all(getattr(obj, property) == value)
        else:
            assert np.allclose(getattr(obj, property), value, atol=0)


def test_staging(datamanager, variant, nopt):
    """Test staging and unstaging."""
    create_input = getattr(datamanager, "create_" + variant)
    obj = create_input()
    properties = obj.export()
    obj.stage()
    for property, value in properties.items():
        if hasattr(value, "units"):
            assert np.allclose(
                getattr(obj, property), value, atol=0*value.units
                )
        elif isinstance(value, dict):
            assert getattr(obj, property) == value
        elif isinstance(value[0], str):
            assert all(getattr(obj, property) == value)
        else:
            assert np.allclose(getattr(obj, property), value, atol=0)
    obj.unstage()
