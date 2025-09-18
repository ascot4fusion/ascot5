import sys
import types
from datetime import datetime

import unyt
import pytest
import numpy as np

from a5py import utils


def test_format2universaldate_correct_format():
    """The function should format datetime into ASCOT5's universal format.

    This is "YYYY-MM-DD HH:MM:SS" with zero padding where appropriate.
    """
    dt = datetime(2023, 5, 17, 8, 3, 45)
    result = utils.format2universaldate(dt)
    assert result == "2023-05-17 08:03:45"


def test_validate_variables_multiple_passes():
    """Multiple variables, shapes, and dtypes should be validated together.

    The return must be a tuple of cast variables in correct order."""
    a = np.array([1, 2])
    b = np.array([[3.0, 4.0]])
    validated = utils.validate_variables(
        [a, b],
        ["a", "b"],
        [(2,), (1, 2)],
        ["int64", "float64"],
    )
    assert isinstance(validated, tuple)
    assert validated[0].shape == (2,)
    assert validated[0].dtype == np.int64
    assert validated[1].shape == (1, 2)
    assert validated[1].dtype == np.float64


def test_incorrect_shape_raises():
    """If a variable cannot be reshaped into the expected shape, the function
    must raise ValueError with a descriptive message.
    """
    arr = np.array([1, 2, 3])
    with pytest.raises(ValueError, match="incompatible shape"):
        utils.validate_variables(arr, "x", (2, 2), "int64")


def test_incorrect_dtype_raises():
    """If a variable cannot be safely cast to the expected dtype without
    changing values, the function must raise ValueError."""
    arr = np.array([1.5, 2.5])  # floats cannot be cast to int64 safely
    with pytest.raises(ValueError, match="incompatible type"):
        utils.validate_variables(arr, "x", (2,), "int64")


def test_scalar_variable_casts_correctly():
    """Scalars (0-dim arrays) must still be validated against shape and dtype.

    They should round-trip correctly when the dtype matches."""
    arr = np.array(42)
    validated = utils.validate_variables(arr, "val", (), "int64")
    assert isinstance(validated, list)
    assert validated[0].shape == ()
    assert validated[0].dtype == np.int64


def test_decorate():
    """Test decoration."""
    s = "hello"
    result = utils.decorate(s)
    assert utils.undecorate(result) == s
    assert result.endswith("\033[0m")


def test_undecorate():
    """Test undecoration."""
    decorated = "\033[01m\033[32mHello\033[0mWorld\033[04m!\033[0m"
    assert utils.undecorate(decorated) == "HelloWorld!"

def test_combination_of_effects():
    """Test that all decorations are in place."""
    s = "combo"
    result = utils.decorate(s, color="green", bold=True, underline=True)
    assert "\033[32m" in result
    assert "\033[01m" in result
    assert "\033[04m" in result
    assert result.endswith("\033[0m")
    assert utils.undecorate(result) == s


def test_too_few_points_raises():
    """A valid abscissa must have at least two points.

    If only one is given, the function should raise ValueError immediately.
    """
    with pytest.raises(ValueError, match="must have at least two points"):
        utils.check_abscissa(np.array([1.0]) * unyt.m, "x")


def test_not_increasing_raises():
    """The abscissa must be strictly increasing.

    If a decrease is present, the function should raise ValueError.
    """
    abscissa = np.array([0.0, 1.0, 0.5]) * unyt.m
    with pytest.raises(ValueError, match="must be increasing"):
        utils.check_abscissa(abscissa, "x")


def test_uniform_passes_and_nonuniform_raises():
    """When uniform=True, differences between consecutive values must match
    the expected uniform spacing.

    Non-uniform abscissas should fail.
    """
    abscissa_uniform = np.linspace(0, 10, 11) * unyt.m
    utils.check_abscissa(abscissa_uniform, "x", uniform=True)

    abscissa_nonuniform = np.array([0, 1, 3, 6, 10]) * unyt.m
    with pytest.raises(ValueError, match="must be uniform"):
        utils.check_abscissa(abscissa_nonuniform, "x", uniform=True)


def test_periodic_passes_and_nonperiodic_raises():
    """When periodic=True, the domain must fit an integer number of 360°.

    A properly periodic grid should pass, while a mismatched domain fails.
    """
    abscissa = np.linspace(0, 360, 5)[:-1] * unyt.deg
    utils.check_abscissa(abscissa, "theta", periodic=True)

    # non-periodic grid: slightly offset
    abscissa_bad = np.linspace(0, 360, 5) * unyt.deg
    with pytest.raises(ValueError, match="must be periodic"):
        utils.check_abscissa(abscissa_bad, "theta", periodic=True)


def test_existing_module():
    """Test that any available module is imported."""
    math = utils.OptionalDependency("math")
    assert math.sqrt(4) == 2


def test_missing_module():
    """Test that *using* (not importing) a missing module raises exception."""
    fake_module = utils.OptionalDependency("fake_module")
    with pytest.raises(ImportError):
        _ = fake_module.some_function()


def test_lazy_loading(monkeypatch):
    """Test that the module is only imported on usage and only on the first
    time.
    """
    called = {}
    counter = {"count": 0}
    def fake_import(name):
        called["name"] = name
        counter["count"] += 1
        return types.SimpleNamespace(foo=lambda: "bar")

    monkeypatch.setitem(sys.modules, "arnoldsbigpackage", None)
    monkeypatch.setattr("importlib.import_module", fake_import)

    opt = utils.OptionalDependency("arnoldsbigpackage")
    assert called == {} # Not called yet
    opt.foo()
    assert called["name"] == "arnoldsbigpackage" # After the first call
    opt.foo()
    assert counter["count"] == 1 # Imported only once


@pytest.fixture()
def node():
    """Create an empty ImmutableNode."""
    return utils.ImmutableStorage()


def test_initialization(node):
    """Test that the initialized node is unfrozen and empty."""
    assert not node._frozen
    assert not len(node._tags)
    assert node._active is None


def test_freeze_and_unfreeze(node):
    """Test that the node is immutable when frozen and mutable when unfrozen."""
    VALUE = object()
    node._freeze()
    assert node._frozen
    with pytest.raises(AttributeError):
        node.new_attr = VALUE
    with pytest.raises(AttributeError):
        node.new_attr

    node._unfreeze()
    assert not node._frozen
    node.new_attr = VALUE
    assert node.new_attr is VALUE


def test_modify_attributes_context(node):
    """Test that the node can can be modified inside modify_attributes context.
    """
    VALUE = object()
    node._freeze()
    with node._modify_attributes():
        node.new_attr = VALUE
    assert node.new_attr is VALUE


def test_setgetitem(node):
    """Test that the node attributes can be assigned and referenced both with
    getattr and getitem methods, and set with either setattr or setitem
    (unless the node is frozen in which case neither should work).
    """
    VALUE, ANOTHER_VALUE, YET_ANOTHER_VALUE = object(), object(), object()
    node["attr"] = VALUE
    assert node.attr is VALUE
    assert node["attr"] is VALUE

    node.attr = ANOTHER_VALUE
    assert node.attr is ANOTHER_VALUE
    assert node["attr"] is ANOTHER_VALUE

    node._freeze()
    with pytest.raises(AttributeError):
        node.attr = YET_ANOTHER_VALUE
    with pytest.raises(AttributeError):
        node["attr"] = YET_ANOTHER_VALUE
    assert node.attr is ANOTHER_VALUE
    assert node["attr"] is ANOTHER_VALUE
