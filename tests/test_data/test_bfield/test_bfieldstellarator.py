"""Test bfieldstellarator.py."""

import pytest
import numpy as np
import unyt

from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data(request, ascot):
    """Test data."""


@pytest.mark.skip("Not yet implemented.")
@pytest.mark.parametrize(
    "ascot",
    [
        pytest.param(True, id="with-file"),
        pytest.param(False, id="no-file"),
    ],
    indirect=True,
)
@pytest.mark.parametrize(
    "data, optional",
    [
        pytest.param(False, False, id="required-data"),
        pytest.param(True, True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldstellarator_init(ascot, data, optional):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_eval_axis_rz(ascot, data, helpers):
    """Test that magnetic axis evaluationr returns correct values."""


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_eval_psi_rho_and_gradients(ascot, data, helpers):
    """Test values and gradients of psi and rho."""


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_value_of_b(ascot, data, helpers):
    """Test that the value of magnetic field is correct."""


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_eval_b_and_bdb(ascot, data, helpers):
    """Test evaluation of magnetic field derivatives."""


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_eval_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""


@pytest.mark.skip("Not yet implemented.")
def test_bfieldstellarator_periodicity(ascot, data, helpers):
    """Test that the field is periodic when interpolated."""
