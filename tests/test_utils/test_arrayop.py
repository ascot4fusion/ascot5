"""Tests for a5py.utils.arrayop."""
import unyt
import pytest
import numpy as np

from a5py import utils

@pytest.mark.parametrize(
    "arr", [
        np.full((3, 2), 1.0)*unyt.m,
        np.full((2, 3), 1.0)*unyt.m,
        np.full((3, 2), 0.001)*unyt.km,
        np.full((2, 3), 1, dtype="i4")*unyt.m,
        None,
    ],
    ids=["as-is", "transpose", "compatible-units", "convert-type", "default",]
    )
def test_validate_variables(arr):
    """Test that validate_variables() works as intended when compatible
    variables are passed.
    """
    default = np.full((3, 2), 1.0) * unyt.m
    with utils.validate_variables() as v:
        validated = v.validate(
            "x", arr, (3, 2), "m", default=default.v,
            )
    assert isinstance(validated, unyt.unyt_array)
    assert validated.shape == (3, 2)
    assert validated.dtype == np.float64
    assert str(validated.units) == "m"
    arr = default if arr is None else arr
    assert np.isclose(validated.ravel(), arr.ravel()).all()


@pytest.mark.parametrize(
    "arr, err",[
        (np.full((2,), 1)*unyt.m, "incompatible shape"),
        ("peruna ", "must be int, float, list"),
        (np.full((2, 2), 1)*unyt.e, "incompatible units"),
        (None, "Default value must be provided"),
        ],
    ids=["shape", "dtype", "units", "no-default",]
    )
def test_validate_variables_invalid(arr, err):
    """Test that using validate_variables() with invalid variables raises
    ValueError.
    """
    with pytest.raises(ValueError, match=err):
        with utils.validate_variables() as v:
            v.validate("x", arr, (2, 2), dtype="f8", units="m")


def test_validate_variables_no_units():
    """Test that validate_variables() issues warning when argument is assumed
    to be in correct units when no units are passed except when the assumed
    units are dimensionless.
    """
    with pytest.warns(UserWarning, match="given without dimensions"):
        with utils.validate_variables() as v:
            arr = v.validate("x", np.full((2, 2), 1), (2, 2), "m")
    assert str(arr.units) == "m"
    with utils.validate_variables() as v:
        arr = v.validate("x", np.full((2, 2), 1), (2, 2), "1")
    assert str(arr.units) == "dimensionless"


@pytest.mark.parametrize(
    "shape, expected", [
        ((), (),),
        ((1,), (),),
        ((1, 1), (),),
        ((5,), (-1,),),
        ((5, 1), (-1,),),
        ],
)
def test_validate_variables_1d_arrays(shape, expected):
    """Test that validate_variables() works with special 1D array cases."""
    with utils.validate_variables() as v:
        arr = v.validate("x", np.full(shape, 1), expected)
    if expected == (-1,):
        expected = (np.multiply.reduce(shape),)
    assert arr.shape == expected


def test_validate_variables_invalid_1d():
    """Test that "shape=(-1,)" option only accepts 1D arrays."""
    with pytest.raises(ValueError, match="expected 1D"):
        with utils.validate_variables() as v:
            v.validate("x", np.full((2, 2), 1), (-1,))


def test_size():
    """Test that size() works as intended."""
    assert utils.size(np.array([1, 2, 3])) == 3
    assert utils.size(1*unyt.m) == 1
    assert utils.size([]) == 0


def test_scalar2array():
    """Test that scalar2array() works as intended."""
    assert np.isclose(utils.scalar2array(1, (3, 2)), np.full((3, 2), 1)).all()
    assert np.isclose(utils.scalar2array(1*unyt.m, (3, 2)), np.full((3, 2), 1)*unyt.m).all()
    assert utils.scalar2array(None, (3, 2)) is None
    assert np.isclose(utils.scalar2array(np.full((3, 2), 1), (3,2)), np.full((3, 2), 1)).all()


def test_too_few_points_raises():
    """A valid abscissa must have at least two points.

    If only one is given, the function should raise ValueError immediately.
    """
    with pytest.raises(ValueError, match="must have at least two points"):
        utils.validate_abscissa(np.array([1.0]) * unyt.m, "x")


def test_not_increasing_raises():
    """The abscissa must be strictly increasing.

    If a decrease is present, the function should raise ValueError.
    """
    abscissa = np.array([0.0, 1.0, 0.5]) * unyt.m
    with pytest.raises(ValueError, match="must be increasing"):
        utils.validate_abscissa(abscissa, "x")


def test_uniform_passes_and_nonuniform_raises():
    """When uniform=True, differences between consecutive values must match
    the expected uniform spacing.

    Non-uniform abscissas should fail.
    """
    abscissa_uniform = np.linspace(0, 10, 11) * unyt.m
    utils.validate_abscissa(abscissa_uniform, "x", uniform=True)

    abscissa_nonuniform = np.array([0, 1, 3, 6, 10]) * unyt.m
    with pytest.raises(ValueError, match="must be uniform"):
        utils.validate_abscissa(abscissa_nonuniform, "x", uniform=True)


def test_periodic_passes_and_nonperiodic_raises():
    """When periodic=True, the domain must fit an integer number of 360°.

    A properly periodic grid should pass, while a mismatched domain fails.
    """
    abscissa = np.linspace(0, 360, 5)[:-1] * unyt.deg
    utils.validate_abscissa(abscissa, "theta", periodic=True)

    # non-periodic grid: slightly offset
    abscissa_bad = np.linspace(0, 360, 5) * unyt.deg
    with pytest.raises(ValueError, match="must be periodic"):
        utils.validate_abscissa(abscissa_bad, "theta", periodic=True)
