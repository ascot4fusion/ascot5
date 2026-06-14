import pytest
import unyt
import numpy as np

from a5py.data.hist import Dist


@pytest.fixture
def dist():
    bin_edges = {"r": unyt.unyt_array(np.linspace(0, 1, 11), unyt.m)}
    dist = Dist(np.full((10, 11), 0.0), bin_edges)
    return dist


@pytest.mark.parametrize(
    "bin_edges, values",
    [
        pytest.param(
            {
                "x1": unyt.unyt_array(np.linspace(0, 1, 11), unyt.m),
                "x2": unyt.unyt_array(np.linspace(0, 1, 12), unyt.s),
            },
            np.full((10, 11), 1.0),
            id="vanilla",
        ),
        pytest.param(
            {"x1": unyt.unyt_array(np.linspace(0, 1, 11), unyt.m)},
            np.full((10,), 1.0),
            id="one-dimensional",
        ),
        pytest.param(
            {"x1": np.linspace(0, 1, 11)}, np.full((10,), 1.0), id="dimensionless-axis"
        ),
    ],
)
def test_dist_init(bin_edges, values):
    dist = Dist(values, bin_edges)
    assert np.isclose(dist.histogram, values).all()
    for k in bin_edges:
        assert np.isclose(dist.bin_edges[k], bin_edges[k]).all()

    vol = unyt.unyt_array([1.0])
    units = unyt.unyt_array([1.0])
    for k in dist.bin_edges:
        units /= dist.bin_edges[k].units
        vol *= np.diff(dist.bin_edges[k][:2])

    assert dist.distribution.units.same_dimensions_as(units.units)
    assert dist.distribution.ravel()[0] == 1.0 / vol


def test_dist_init_with_species():
    dist = Dist(
        histogram=np.full((10,), 0.0),
        bin_edges={"x": np.linspace(0, 1, 11)},
        species="H",
        charge=1,
    )
    assert dist.species.name == "H"
    assert dist.charge == 1


@pytest.mark.parametrize(
    "bin_edges, values, message",
    [
        pytest.param(
            {"x": np.full((11, 11), 0.0)},
            np.full((10,), 0.0),
            "must be 1-dimensional",
            id="axis-not-1d",
        ),
        pytest.param(
            {"x": np.array([0])},
            np.full((1,), 0.0),
            "must contain at least two entries",
            id="axis-size-not-geq-2",
        ),
        pytest.param(
            {"x": np.array([0, 1, 3])},
            np.full((10,), 0.0),
            "is not uniformly spaced",
            id="axis-not-uniform",
        ),
        pytest.param(
            {"x": np.linspace(0, 1, 4)},
            np.full((2,), 0.0),
            "must have length 3",
            id="axis-size-not-same-as-in-values",
        ),
        pytest.param(
            {"x": np.linspace(0, 1, 11)},
            np.full((10, 5), 0.0),
            "edges must have the same number of dimensions as the histogram",
            id="number-of-dimensions-not-same-as-in-values",
        ),
    ],
)
def test_dist_init_raises(bin_edges, values, message):

    with pytest.raises(ValueError, match=message):
        Dist(values, bin_edges)


def test_dist_init_species_requires_charge():
    with pytest.raises(ValueError, match="`charge` must be specified"):
        Dist(np.full((10,), 0.0), {"x": np.linspace(0, 1, 11)}, species="H")


def test_dist_copy():
    dist = Dist(np.full((10,), 0.0), {"x": np.linspace(0, 1, 11)}, "H", charge=1)
    dist._distribution[0] = 1.0
    copy = dist.copy()
    assert np.isclose(copy.distribution, dist.distribution).all()
    for k in dist.bin_edges:
        assert np.isclose(copy.bin_edges[k], dist.bin_edges[k]).all()

    assert dist.distribution is not copy.distribution
    dist._distribution[0] = 2.0
    assert copy.distribution[0] == 1.0

    assert dist.species == copy.species
    assert dist.charge == copy.charge


@pytest.mark.parametrize(
    "multiplier, dimensions, result",
    [
        pytest.param(
            np.linspace(0, 1, 1), (0,), np.full((1, 2, 3, 2), 0.0), id="first-axis"
        ),
        pytest.param(
            np.linspace(0, 1, 2),
            (1,),
            np.array(
                [
                    [
                        [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                        [[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]],
                    ]
                ]
            ),
            id="second-axis",
        ),
        pytest.param(
            np.linspace(0, 1, 3)*unyt.m,
            (2,),
            np.array(
                [
                    [
                        [[0.0, 0.0], [1.5, 2.0], [5.0, 6.0]],
                        [[0.0, 0.0], [4.5, 5.0], [11.0, 12.0]],
                    ]
                ]
            )*unyt.m,
            id="third-axis",
        ),
        pytest.param(
            2.0*unyt.m,
            None,
            np.array(
                [
                    [
                        [[2.0, 4.0], [6.0, 8.0], [10.0, 12.0]],
                        [[14.0, 16.0], [18.0, 20.0], [22.0, 24.0]],
                    ]
                ]
            )*unyt.m,
            id="scalar",
        ),
        pytest.param(
            np.array(
                [
                    [
                        [[12.0, 11.0], [10.0, 9.0], [8.0, 7.0]],
                        [[6.0, 5.0], [4.0, 3.0], [2.0, 1.0]],
                    ]
                ]
            ),
            None,
            np.array(
                [
                    [
                        [[12., 22.], [30., 36.], [40., 42.]],
                        [[42., 40.], [36., 30.], [22., 12.]],
                    ]
                ]
            ),
            id="full-matrix",
        ),
        pytest.param(
            np.array([[1.0, 2.0, 1.0], [2.0, 1.0, 2.0]]),
            (1, 2),
            np.array([[[1.0, 4.0, 3.0], [8.0, 5.0, 12.0]]]),
            id="partial-matrix",
        ),
    ],
)
def test_dist_multiply(multiplier, dimensions, result):
    bin_edges = {
        "x1": np.linspace(0, 1, 2),
        "x2": np.linspace(0, 1, 3),
        "x3": np.linspace(0, 1, 4),
        "x4": np.linspace(0, 1, 3),
    }

    shape = []
    for dim in bin_edges.values():
        shape.append(dim.size - 1)

    histogram = np.ones(shape)
    histogram += np.arange(histogram.size).reshape(histogram.shape)
    dist = Dist(histogram, bin_edges)
    dist.multiply(multiplier, dimensions)
    print(dist.histogram)
    assert np.isclose(dist.histogram, result).all()

def test_dist_multiply_raises():
    pass


@pytest.mark.parametrize(
    ("slices", "expected_edges"),
    [
        (
            {"x1": 1},
            {"x1": np.arange(4)[1:3], "x2": np.arange(5), "x3": np.arange(6)},
        ),
        (
            {"x1": -1},
            {"x1": np.arange(4)[-2:], "x2": np.arange(5), "x3": np.arange(6)},
        ),
        (
            {"x1": slice(1, 3)},
            {"x1": np.arange(4)[1:4], "x2": np.arange(5), "x3": np.arange(6)},
        ),
        (
            {"x1": 1, "x3": -1},
            {"x1": np.arange(4)[1:3], "x2": np.arange(5), "x3": np.arange(6)[-2:]},
        ),
    ],
    ids=["slice-integer", "slice-end", "slice-slice", "two-slices"],
)
def test_slice_valid(slices, expected_edges):
    bin_edges = {
        "x1": np.arange(4),
        "x2": np.arange(5),
        "x3": np.arange(6),
    }

    # Construct histograms from bin edges
    histogram = np.einsum(
        "i,j,k->ijk",
        bin_edges["x1"][1:],
        bin_edges["x2"][1:]*10,
        bin_edges["x3"][1:]*100,
        )
    expected_histogram = np.einsum(
        "i,j,k->ijk",
        expected_edges["x1"][1:],
        expected_edges["x2"][1:]*10,
        expected_edges["x3"][1:]*100,
        )
    dist = Dist(histogram, bin_edges)
    dist.slice(slices)

    assert dist.distribution.shape == expected_histogram.shape
    assert np.isclose(dist.distribution, expected_histogram).all()
    for dimension, edges in expected_edges.items():
        assert dimension in dist.bin_edges
        assert dist.bin_edges[dimension].shape == edges.shape
        assert np.isclose(dist.bin_edges[dimension], edges).all()

@pytest.mark.parametrize(
    ("slices", "error", "message"),
    [
        ({"x1": 4}, ValueError, "index 4 is out of bounds for axis 0 with size 3"),
        ({"x1": -4}, ValueError, "index -4 is out of bounds for axis 0 with size 3"),
        ({"x4": 1}, KeyError, "Unknown dimension(s): 'x4'"),
        ({"x1": (3,4)}, TypeError, "Slice for dimension 'x1' must be int or slice"),
    ]
)
def test_dist_slice_raises(slices, error, message):
    dist = Dist(
        np.full((3,4,5), 0),
        {
            "x1": np.arange(4),
            "x2": np.arange(5),
            "x3": np.arange(6),
        }
    )
    with pytest.raises(error, match=message):
        dist.slice(slices)


def test_dist_integrate():
    dist = Dist(
        np.full((3,4,5), 1),
        {
            "x1": np.arange(4),
            "x2": np.arange(5),
            "x3": np.arange(6),
        }
    )
    dist.integrate({"x1": None, "x2": None, "x3": None})
    print(dist.distribution)

def test_dist_integrate_raises():
    pass

@pytest.mark.parametrize(
    "coords",
    [
        {"x": 1.5, "y": 0.},
        {"x": 1.0, "y": 1.0},
        {"x": 0., "y": 2.},
        {"x": 0.1,},
        {"y": 0.1,},
    ],
)
def test_interpolate(coords):
    x = np.linspace(0., 3., 4)
    y = np.linspace(0., 4., 5)
    xc, yc = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:])

    # f(x,y) = x + 2y
    xx, yy = np.meshgrid(xc, yc, indexing="ij")
    data = xx + 2*yy

    dist = Dist(
        histogram=data,
        bin_edges={
            "x": x,
            "y": y,
        },
    )
    result = dist.interpolate(coords)

    xq = xc if not "x" in coords else coords["x"]
    yq = yc if not "y" in coords else coords["y"]
    expected = np.asarray(xq + 2*yq)

    assert result.shape == expected.shape
    assert np.isclose(result, expected).all()

@pytest.mark.parametrize(
    "coords",
    [
        {"x": 1.5, "y": 0.},
        {"x": np.array([0.5, 1.5]), "y": np.array([1.0, 1.0])},
    ],
)
def test_interpolate_points(coords):
    x = np.linspace(0., 3., 4)
    y = np.linspace(0., 4., 5)
    xc, yc = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:])

    # f(x,y) = x + 2y
    xx, yy = np.meshgrid(xc, yc, indexing="ij")
    data = xx + 2*yy

    dist = Dist(
        histogram=data,
        bin_edges={
            "x": x,
            "y": y,
        },
    )
    xq = np.array([coords["x"], coords["y"]]).T
    result = dist.interpolate_points(xq)

    xq = xc if not "x" in coords else coords["x"]
    yq = yc if not "y" in coords else coords["y"]
    expected = np.asarray(xq + 2*yq)

    assert result.shape == expected.shape
    assert np.isclose(result, expected).all()

def test_dist_interpolate_raises():
    pass
