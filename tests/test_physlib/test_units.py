import pytest
import unyt
import numpy as np

import a5py.physlib as physlib
from a5py.exceptions import AscotUnitWarning

@pytest.fixture(
        name="value",
        params=[
            "float", "int", "[float]", "[int]", "ndarray([int])",
            "ndarray([float])", "unyt_array([int])", "unyt_array([float])",
            ],
        )
def fixture_value(request):
    """Fixture providing arrays of different types."""
    match request.param:
        case "float":
            return 1.0
        case "int":
            return 1
        case "[float]":
            return [1.0]
        case "[int]":
            return [1]
        case "ndarray([int])":
            return np.array([1, 2], dtype=int)
        case "ndarray([float])":
            return np.array([1.0, 2.0], dtype=float)
        case "unyt_array([int])":
            return 1 * unyt.inch
        case "unyt_array([float])":
            return 1.0 * unyt.inch


def test_match_units(value):
    """Test match_units function with default parameters."""
    out = physlib.match_units(value, unit="m", strip=False, auto_assign=True)
    assert value is not out
    assert out.units == unyt.m


def test_match_units_without_auto_assign(value):
    """Test match_units function with auto_assign=False."""
    if hasattr(value, "units"):
        out = physlib.match_units(value, unit="m", strip=False, auto_assign=False)
        assert value is not out
        assert out.units == unyt.m
    else:
        with pytest.raises(ValueError):
            physlib.match_units(value, unit="m", strip=False, auto_assign=False)


def test_match_units_with_strip(value):
    """Test match_units function with strip=False."""
    out = physlib.match_units(value, unit="m", strip=True, auto_assign=True)
    assert value is not out
    assert not hasattr(out, "units")


def test_parse_units():
    """Test the parse_units decorator."""
    @physlib.parseunits(arg_units="m", kwarg_units="m")
    def func(
        arg_units, arg_not_listed, arg_units_not_listed,
        kwarg_units=None, kwarg_not_listed=None, kwarg_units_not_listed=None,
        ):
        assert arg_units.units == unyt.m
        assert not hasattr(arg_not_listed, "units")
        assert arg_units_not_listed.units == unyt.inch
        assert kwarg_units.units == unyt.m
        assert not hasattr(kwarg_not_listed, "units")
        assert kwarg_units_not_listed.units == unyt.inch
    arg1, arg2, arg3 = 1*unyt.inch, 1, 1*unyt.inch
    kwarg1, kwarg2, kwarg3 = 1*unyt.inch, 1, 1*unyt.inch
    func(arg1, arg2, arg3, kwarg_units=kwarg1, kwarg_not_listed=kwarg2,
         kwarg_units_not_listed=kwarg3)
    assert arg1.units == unyt.inch
    assert not hasattr(arg2, "units")
    assert arg3.units == unyt.inch
    assert kwarg1.units == unyt.inch
    assert not hasattr(kwarg2, "units")
    assert kwarg3.units == unyt.inch


def test_parse_units_with_strip():
    """Test the parse_units decorator with strip=True."""
    @physlib.parseunits(arg="m", kwarg="m", strip=True)
    def func(arg, kwarg=None):
        assert not hasattr(arg, "units")
        assert not hasattr(kwarg, "units")
    arg, kwarg = 1*unyt.inch, 1*unyt.inch
    func(arg, kwarg=kwarg)
    assert arg.units == unyt.inch
    assert kwarg.units == unyt.inch


def test_parse_units_incompatible_dimensions():
    """Test the parse_units decorator when argument has wrong dimensionality."""
    @physlib.parseunits(arg="m", kwarg="m")
    def func(arg, kwarg=None):
        pass
    with pytest.raises(ValueError):
        func(1*unyt.s, kwarg=1*unyt.m)
    with pytest.raises(ValueError):
        func(1*unyt.m, kwarg=1*unyt.s)


def test_parse_units_assumed_units():
    """Test the parse_units decorator when units are assumed."""
    @physlib.parseunits(arg="m", kwarg="m")
    def func(arg, kwarg=None):
        pass
    with pytest.warns(AscotUnitWarning):
        func(1, kwarg=1*unyt.m)
    with pytest.warns(AscotUnitWarning):
        func(1*unyt.m, kwarg=1)
