"""Test options."""
import os

import pytest
import numpy as np

from a5py.data import AscotData
from a5py.data import options
from a5py.data.options import constraints


@pytest.fixture(name="datamanager")
def fixture_datamanager():
    yield AscotData(hdf5file=("testascot.h5", False))
    os.unlink("testascot.h5")


@pytest.mark.parametrize(
    "store_hdf5",
    [True, False],
    ids=["store_hdf5", "dont_store_hdf5"],
    )


def single_line():
    """Single line: (> 0), default=1."""

def two_line():
    """Multi-line:
    (> 0), default=1.
    """

def multi_line():
    """Multi-line split at
    constraint:
    (> 0),
    default=1.
    """

@pytest.mark.parametrize(
    "fun",
    [single_line, two_line, multi_line],
    ids=["single_line", "two_line", "multi_line"],
)
def test_summarize_property(fun):
    info = constraints.summarize_property(fun)
    assert info["name"] == fun.__name__
    assert info["constraint"] == "(> 0)"
    assert info["default"] == "1"


@pytest.mark.parametrize(
    "constraint, valid_values, invalid_values",
    [
        ("{1, 2, 3, 4}", (1, 2, 3, 4), (0, 5, 1.5)),
        ("(x)", (1,), ()),
        ("(> 0)", (1,), (0,)),
        ("[a, b]", ([1, 2], (1, 2)), ([], [1, 2, 3])),
        ("[a > 0, b > 0]", ([1, 1],), ([0, 1], [1, 0], [0, 0])),
        ("[0 <= a < 360, 0 < b <= 360]", ([1, 2],), ([-1, 2], [1, 0], [360, 1], [1, 361])),
        ("[0 <= a < 360, ...]", (0, [0], [1, 359],), ([0, 361],)),
    ],
    ids=[
        "sequence", "no-constraint", "positive", "2-array-no-limits",
        "2-array-lower-limit", "2-array-both-limits", "unlimited-array-with-limits",
    ]
    )
def test_parse(constraint, valid_values, invalid_values):
    for value in valid_values:
        assert constraints.enforce_constraint(value, constraint)

    for value in invalid_values:
        assert not constraints.enforce_constraint(value, constraint)


@pytest.mark.parametrize(
    "store_hdf5",
    [True, False],
    ids=["store_hdf5", "dont_store_hdf5"],
    )
def test_create(datamanager, store_hdf5):
    """Test input creation."""
    obj = datamanager.create_options(dryrun=True)
    parameters = obj.export()

    obj = datamanager.create_options(**parameters, store_hdf5=store_hdf5)
    assert obj.simulation.simulation_mode == parameters["simulation_mode"]


def test_tostring(datamanager):
    obj = datamanager.create_options()
    #print(obj.export_as_string())
    assert True


def test_stage(datamanager):
    obj = datamanager.create_options()
    obj.stage()
    assert False
