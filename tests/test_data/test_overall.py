"""An overall test to verify all input variants can be created, staged, saved
and data accessed.
"""
import os

import unyt
import pytest
import numpy as np

from a5py.physlib import Species
from a5py.data import AscotData

from .conf import test_data, create_input

@pytest.mark.parametrize("variant", list(test_data.keys()))
@pytest.mark.parametrize("optional", [True, False])
@pytest.mark.parametrize("save", [True, False])
def test_create(variant, optional, save):
    """Test that the factory method produces an object whose properties return
    correct values.
    """
    parameters = test_data[variant]
    try:
        if save:
            data = AscotData(("test.h5", False))
        else:
            data = AscotData()

        factorymethod = getattr(data, f"create_{variant.lower()}")

        if optional:
            parameters = {
                param[1:] if param.startswith("_") else param: value
                for param, value in parameters.items()
                }
        else:
            parameters = {
                param: value for param, value in parameters.items()
                if not param.startswith("_")
                }

        obj = factorymethod(**parameters)
        for param, value in parameters.items():
            if param == "labels":
                obj_dict = getattr(obj, param)
                for label, flag in value.items():
                    assert label in obj_dict
                    assert flag == obj_dict[label]
                for label, flag in obj_dict.items():
                    assert label in value
                    assert flag == value[label]
            elif param == "species":
                species = obj.species
                if not isinstance(value, list):
                    value = [value]
                    species = [species]
                for s1, s2 in zip(species, value):
                    assert s1.name == Species.from_string(s2).name
            else:
                assert np.isclose(getattr(obj, param), value).all()
    finally:
        if save:
            os.unlink("test.h5")


@pytest.mark.parametrize("variant", list(test_data.keys()))
def test_create_invalid_shape_or_units(variant):
    """Test that invalid shape or units raise an error.

    Only checks parameters that are numpy or unit arrays.
    """
    parameters = test_data[variant]
    parameters = {
        param[1:] if param.startswith("_") else param: value
        for param, value in parameters.items()
        }
    for name, val in parameters.items():
        if not isinstance(val, np.ndarray) or not isinstance(val, unyt.unyt_array):
            continue

        parameters[name] = val.T
        create_input(variant, **parameters)
        assert True
        parameters[name] = val

        if val.size == 1:
            parameters[name] = np.full(2, val)
        else:
            parameters[name] = val[1:]

        with pytest.raises(ValueError):
            create_input(variant, **parameters)
        parameters[name] = val

        if isinstance(val, unyt.unyt_array):
            parameters[name] = val * unyt.yomama
            with pytest.raises(ValueError):
                create_input(variant, **parameters)

            parameters[name] = val


def test_export():
    """Test that the data can be exported."""


def test_staging():
    """"""


def test_interpolate():
    """"""


def test_options():
    """"""
    data = AscotData()
    data.create_options().stage()
