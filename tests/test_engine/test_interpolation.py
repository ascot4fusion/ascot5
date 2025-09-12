import pytest
import unyt
import numpy as np

from a5py import Ascot
from a5py.templates import PremadeMagneticField, FlatPlasma

from a5py.engine import interpolate


def generate_inputs(inputs):
    """Generate interpolatable inputs data."""
    a5 = Ascot()
    if "bfield" in inputs:
        template = PremadeMagneticField(a5, field="iter-baseline")
        template.create_input()
    if "efield" in inputs:
        a5.data.create_efieldcartesian(exyz=np.array([1., 2., 3.]))
    if "plasma" in inputs:
        template = FlatPlasma(a5)
        template.create_input()
    if "neutral" in inputs:
        pass
    if "boozer" in inputs:
        pass
    if "mhd" in inputs:
        pass
    return a5


@pytest.fixture
def ascot_with_inputs(request):
    """Generate Ascot object with specified inputs."""
    return generate_inputs(request.param)


@pytest.mark.parametrize(
    "units", [None, "m"], ids=["no-units", "with-units"],
)
def test_init_return_array(units):
    """Test that the init_return_array returns an array with correct shape,
    where all elements are initialized to NaN.
    """
    arr = interpolate.init_return_array((2, 3), units=units)
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (2, 3)
    assert np.isnan(arr).all()


@pytest.mark.parametrize("inputs, invalid_shape", [
    ((1*unyt.m, 5*unyt.deg, 2*unyt.m, 3*unyt.s), False),
    ((unyt.unyt_array([1, 2], "m"), 5*unyt.deg, 2*unyt.m, 3*unyt.s), False),
    ((unyt.unyt_array([1, 2], "m"), 5*unyt.deg, unyt.unyt_array([1, 2, 3], "m"), 3*unyt.s), False),
    ((unyt.unyt_array([[1, 2], [1, 2]], "m"), 5*unyt.deg, 2*unyt.m, 3*unyt.s), True),
    ],
    ids=["all-singleton", "one-array", "arrays-of-different-size", "not-1d-array"],
)
def test_process_query_points_withgrid(inputs, invalid_shape):
    if invalid_shape:
        with pytest.raises(ValueError):
            interpolate.process_query_points(*inputs, grid=True)
        return
    arrshape = tuple(var.size for var in inputs)
    r, phi, z, t = interpolate.process_query_points(*inputs, grid=True)
    assert str(phi.units) == "rad"
    assert all(arr.shape == arrshape for arr in (r, phi, z, t))


@pytest.mark.parametrize("inputs, invalid_shape", [
    ((1*unyt.m, 5*unyt.deg, 2*unyt.m, 3*unyt.s), False),
    ((unyt.unyt_array([1, 2], "m"), 5*unyt.deg, 2*unyt.m, 3*unyt.s), False),
    ((unyt.unyt_array([1, 2], "m"), 5*unyt.deg, unyt.unyt_array([1, 2, 3], "m"), 3*unyt.s), True),
    ((unyt.unyt_array([[1, 2], [1, 2]], "m"), 5*unyt.deg, 2*unyt.m, 3*unyt.s), False),
    ],
    ids=["all-singleton", "one-array", "arrays-of-different-size", "not-1d-array"],
)
def test_process_query_points_withoutgrid(inputs, invalid_shape):
    if invalid_shape:
        with pytest.raises(ValueError):
            interpolate.process_query_points(*inputs, grid=False)
        return
    arrshape = inputs[0].shape # In this test r dictates the shape
    r, phi, z, t = interpolate.process_query_points(*inputs, grid=False)
    assert str(phi.units) == "rad"
    assert all(arr.shape == arrshape for arr in (r, phi, z, t))


@pytest.mark.parametrize(
    "ascot_with_inputs",
    [ ("bfield",), ],
    indirect=True)
@pytest.mark.parametrize( "shape", [(), (5,), (5,5)], )
def test_interpolation_shape(ascot_with_inputs, shape):
    a5 = ascot_with_inputs
    r, phi, z, t = tuple(c*np.ones(shape) for c in (6*unyt.m, 0*unyt.deg, 0*unyt.m, 0*unyt.s))
    br = interpolate.interpolate(r, phi, z, t, 0, "br", bfield=a5.data.bfield.active)["br"]
    assert br.shape == shape
    assert np.unique(br).size == 1


@pytest.mark.parametrize(
    "ascot_with_inputs",
    [ ("bfield",), ],
    indirect=True)
@pytest.mark.parametrize( "qnt", [("br",), ("bphi", "br")] )
def test_interpolation_var(ascot_with_inputs, qnt):
    a5 = ascot_with_inputs
    r, phi, z, t = 6.0*unyt.m, 0.0*unyt.deg, 0.0*unyt.m, 0.0*unyt.s
    values = interpolate.interpolate(r, phi, z, t, 0, *qnt, bfield=a5.data.bfield.active)
    return
    if len(qnt) == 1:
        assert isinstance(values, unyt.unyt_array)
    else:
        assert len(values) == len(qnt)
        for q, v in zip(qnt, values):
            val = interpolate.interpolate(r, phi, z, t, 0, q, bfield=a5.data.bfield.active)
            assert val[q] == v


@pytest.mark.parametrize(
    "inputs, available",
    [
        (("bfield",), ("br", "brdr", "psi", "rho",)),
        (("bfield", "efield",), ("br", "brdr", "psi", "rho", "er",)),
        (("bfield", "plasma",), ("br", "brdr", "psi", "rho", "ne", "te",)),
    ])
def test_interpolation_inputs(inputs, available):
    a5 = generate_inputs(inputs)
    input_variants = {k: getattr(a5.data, k).active for k in inputs}
    r, phi, z, t = 6.0*unyt.m, 0.0*unyt.deg, 0.0*unyt.m, 0.0*unyt.s
    for qnt in ["br", "brdr", "psi", "rho", "er", "ne", "te", "n0", "t0", "theta", "zeta", "alpha", "Phi", "br_mhd", "er_mhd", "phi"]:
        if qnt not in available:
            with pytest.raises(ValueError):
                interpolate.interpolate(r, phi, z, t, 0, qnt, **input_variants)
            continue
        else:
            interpolate.interpolate(r, phi, z, t, 0, qnt, **input_variants)


def test_interpolation_nmode():
    pass


@pytest.mark.parametrize(
    "ascot_with_inputs",
    [ ("bfield",), ],
    indirect=True)
@pytest.mark.parametrize(
    "qnt, unit",
    [ ("br", "T"), ("bphi", "T"), ("bz", "T"), ],
)
def test_interpolation_units(ascot_with_inputs, qnt, unit):
    a5 = ascot_with_inputs
    r, phi, z, t = 6.0*unyt.m, 0.0*unyt.deg, 0.0*unyt.m, 0.0*unyt.s
    val = interpolate.interpolate(r, phi, z, t, 0, qnt, bfield=a5.data.bfield.active)
    assert str(val[qnt].units) == unit


@pytest.mark.parametrize(
    "ascot_with_inputs",
    [ ("bfield",), ],
    indirect=True)
def test_interp(ascot_with_inputs):
    a5 = ascot_with_inputs
    qnt = "jnorm"
    r, phi, z, t = 6.0*unyt.m, 0.0*unyt.deg, 0.0*unyt.m, 0.0*unyt.s
    br = interpolate.evaluate(r, phi, z, t, qnt, bfield=a5.data.bfield.active)
    print(br)
    assert False
