import os
import pytest

from a5py.data.hist import *
from a5py.data.options import SimulationOptions

from a5py.data.access.hdf5io import TreeFileManager

@pytest.fixture()
def file():
    datamanager = TreeFileManager(
        root="root",
        filename="testascot.h5",
        file_exists=False,
        input_categories=[],
        )
    datamanager.set_note(name="hist_1", node="root", note="") # Creates node
    yield datamanager.access_data(name="hist_1", node="root")

    os.unlink("testascot.h5")

@pytest.fixture
def hist():
    opt = SimulationOptions.from_dict(
        histograms=[
            {
                "dimensions": [("r", 0, 1, 10)],
                "charge_interval": (1,1),
            }
        ]
    )
    return Hist.from_params(opt.histograms[0])


def test_hist_init_from_params():
    opt = SimulationOptions.from_dict(
        histograms=[
            {
                "dimensions": [("r", 0, 1, 10)],
                "charge_interval": (1,1),
            }
        ]
    )
    hist = Hist.from_params(opt.histograms[0])
    for i, dimension in enumerate(DIMENSIONS.keys()):
        if dimension == "r":
            assert np.isclose(
                hist.dimensions["r"],
                unyt.unyt_array(np.linspace(0, 1, 11), unyt.m),
                ).all()
            assert hist._cdata.dimensions[i].nbin == 10
        else:
            assert dimension not in hist.dimensions
            assert hist._cdata.dimensions[i].nbin == 0


def test_hist_store_hdf5(hist, file):
    hist._cdata.values[3] = 1.0
    values = hist.values
    hist.save(file)
    
    hist2 = Hist()
    hist2._file = file
    assert all(hist2.values == values)


def test_hist_properties(hist):
    assert len(hist.dimensions) == 1
    assert np.isclose(
        hist.dimensions["r"],
        unyt.unyt_array(np.linspace(0, 1, 11), unyt.m)
        ).all()

    with pytest.raises(AttributeError):
        hist.dimensions = {
            "r": unyt.unyt_array(np.linspace(0, 1, 11), unyt.m)
            }

    hist.dimensions["z"] = unyt.unyt_array(np.linspace(0, 1, 11), unyt.m)
    assert "z" not in hist.dimensions

    hist.dimensions["r"][0] = 2.0
    assert hist.dimensions["r"][0].value == 0.0

    assert np.isclose(
        hist.values,
        np.full((10,), 0.),
    ).all()

    with pytest.raises(AttributeError):
        hist.values = np.full((10,), 1)

    with pytest.raises(ValueError):
        hist.values[0] = 1


def test_hist_extract_dist(hist):
    dist = hist.extract(("r",))
    assert np.isclose(
        dist.bin_edges["r"],
        unyt.unyt_array(np.linspace(0, 1, 11), unyt.m)
        ).all()

    assert np.isclose(
        dist.histogram,
        np.full((10,), 0.)
    ).all()





