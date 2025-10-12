"""Tests for a5py.utils.OptionalDependency."""
import sys
import types
import pytest
from typing import TYPE_CHECKING

from a5py.utils import OptionalDependency

def test_existing_module():
    """Test that any available module is imported."""
    math = OptionalDependency("math")
    if TYPE_CHECKING:
        import math
    assert math.sqrt(4) == 2


def test_missing_module():
    """Test that *using* (not importing) a missing module raises exception."""
    fake_module = OptionalDependency("fake_module")
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

    opt = OptionalDependency("arnoldsbigpackage")
    assert called == {} # Not called yet
    opt.foo()
    assert called["name"] == "arnoldsbigpackage" # After the first call
    opt.foo()
    assert counter["count"] == 1, "Module imported more than once."
