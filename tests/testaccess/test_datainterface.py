"""Test accessing data on disk or in C struct with unified interface."""
import os
import ctypes
import pytest

from collections import namedtuple
from unittest.mock import MagicMock

import unyt
import numpy as np

from a5py import utils
from a5py.physlib import parseunits
from a5py.exceptions import AscotIOException

from a5py.data.access.hdf5 import HDF5Manager
from a5py.data.access.variants import InputVariant
from a5py.data.access.dataholder import Format

from .conftest import QID1, DATE, NOTE, INPUTVAR, FNEMPTY

class Pointcloud(InputVariant):
    """A prototype class for testing purposes.

    All actual data classes inherit both from MetaDataHolder and DataHolder.
    Therefore, this class is a prototype for the data classes. This class
    doesn't have any real functionality beyond testing.
    """

    class Struct(ctypes.Structure):
        """Data stored in C."""
        _pack_ = 1
        _fields_ = [
            ("Npoint", ctypes.c_int32),
            ("volume", ctypes.c_double),
            ("dirvec", ctypes.c_double * 3),
            ("offload_array_size", ctypes.c_size_t),
            ("int_offload_array_size", ctypes.c_size_t),
            ]

    def __init__(self):
        super().__init__(
            qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR,
            struct=Pointcloud.Struct(),
            )

    @property
    def npoint(self):
        """Number of points."""
        if self._format == Format.HDF5:
            return self._read_hdf5("xcoords").size
        return self._struct_.Npoint

    @property
    def volume(self):
        """Volume of the cloud."""
        if self._format == Format.HDF5:
            return self._read_hdf5("volume")
        return self._struct_.volume * unyt.m**3

    @property
    def dirvec(self):
        """Cloud direction."""
        if self._format == Format.HDF5:
            return self._read_hdf5("directionvector")
        return self._from_struct_("dirvec")

    @property
    def x(self):
        """Array of x-coordinates."""
        if self._format == Format.HDF5:
            return self._read_hdf5("xcoords")
        return self._slice_offload_array(0, self.npoint) * unyt.m

    @property
    def y(self):
        """Array of y-coordinates."""
        if self._format == Format.HDF5:
            return self._read_hdf5("ycoords")
        slice_ = (self.npoint, 2*self.npoint)
        return self._slice_offload_array(*slice_) * unyt.m

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already in stored in the file.")
        data = self.export()
        data["directionvector"] = data.pop("dirvec")
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        """Return a dictionary with sufficient data to duplicate this instance.
        """
        data = {
            "volume":self.volume,
            "dirvec":self.dirvec,
            "xcoords":self.x,
            "ycoords":self.y,
        }
        return data

    @classmethod
    def init_ctypes(cls, volume, dirvec, xcoords, ycoords):
        """Initialize point cloud from give nd.array's."""
        obj = cls()
        obj._format = Format.CSTRUCT

        obj._struct_.Npoint = xcoords.size
        obj._struct_.volume = volume
        obj._struct_.dirvec[:] = dirvec

        offload_array_size = xcoords.size + ycoords.size
        obj._offload_array_ = (ctypes.c_double * offload_array_size)()
        obj._struct_.offload_array_size = offload_array_size

        obj._append_offload_array(xcoords, ycoords)
        return obj


@parseunits(volume="m**3", dirvec="1", xcoords="m", ycoords="m", strip=True)
def create_pointcloud(volume, dirvec, xcoords, ycoords):
    """Creates a point cloud instance.

    This mimics the factory methods creating the actual data objects. Inputs
    should be checked in these methods.
    """
    ncoord = xcoords.size
    volume, dirvec, xcoords, ycoords = utils.validate_variables(
        [volume, dirvec, xcoords, ycoords],
        ["volume", "dirvec", "xcoords", "ycoords"],
        [(), (3,), (ncoord,), (ncoord,)],
        "f8"
        )
    return Pointcloud.init_ctypes(volume, dirvec, xcoords, ycoords)


@pytest.fixture
def testdata():
    """Test data for point cloud initialization."""
    Data = namedtuple("data", ["xcoords", "ycoords", "dirvec", "volume"])
    return Data(
        np.array([1.0, 2.0, 3.0, 4.0]),
        np.array([-1.0, -2.0, -3.0, -4.0]),
        np.array([1/np.sqrt(2), -1/np.sqrt(2), 0]),
        -5.0,
        )


@pytest.fixture
def testdatadict(testdata):
    return {
        "volume":testdata.volume*unyt.m**3,
        "dirvec":testdata.dirvec,
        "xcoords":testdata.xcoords*unyt.m,
        "ycoords":testdata.ycoords*unyt.m,
        }


@pytest.fixture
def hdf5manager():
    yield HDF5Manager(FNEMPTY, file_exists=False)
    os.unlink(FNEMPTY)


def test_init_ctypes(testdata):
    """Test initialization in ctypes."""
    pc = Pointcloud.init_ctypes(
        testdata.volume, testdata.dirvec, testdata.xcoords, testdata.ycoords
        )
    assert pc.npoint == 4
    assert pc.volume.v == testdata.volume
    assert np.allclose(pc.x.v, testdata.xcoords)
    assert np.allclose(pc.y.v, testdata.ycoords)
    assert np.allclose(pc.dirvec, testdata.dirvec)

    value_to_ensure_data_is_copied = 10.0
    testdata.dirvec[0] = value_to_ensure_data_is_copied
    assert pc.dirvec[0] != value_to_ensure_data_is_copied
    testdata.xcoords[0] = value_to_ensure_data_is_copied
    assert pc.x[0].v != value_to_ensure_data_is_copied


def test_struct_data_is_immutable(testdata):
    """Test that the getter values are copies and setters are unavailable."""
    value_to_ensure_data_is_copied = 10.0
    pc = Pointcloud.init_ctypes(
        testdata.volume, testdata.dirvec, testdata.xcoords, testdata.ycoords
        )
    pc.dirvec[0] = value_to_ensure_data_is_copied
    assert pc.dirvec[0] != value_to_ensure_data_is_copied
    copied_value = pc.dirvec[0]
    pc._struct_.dirvec[0] = value_to_ensure_data_is_copied
    assert pc.dirvec[0] == value_to_ensure_data_is_copied
    assert copied_value != value_to_ensure_data_is_copied
    copied_value = pc.x[0]
    pc._offload_array_[0] = value_to_ensure_data_is_copied
    assert pc.x[0].v == value_to_ensure_data_is_copied
    assert copied_value != value_to_ensure_data_is_copied

    with pytest.raises(AttributeError):
        pc.volume = np.array([1.0]) * unyt.m**3
    with pytest.raises(AttributeError):
        pc.dirvec = np.array([1,0,0])
    with pytest.raises(AttributeError):
        pc.x = np.array([1.0, 2.0, 3.0, 4.0]) * unyt.m


def test_create(testdatadict):
    """Test creating point cloud from dictionary."""
    testdatadict["volume"] *= 1.0 * unyt.cm**3 / unyt.m**3
    pc = create_pointcloud(**testdatadict)
    assert pc.npoint == 4
    assert np.allclose(pc.volume.v, -5.0e-6)
    assert np.allclose(pc.x.v, testdatadict["xcoords"].v)
    assert np.allclose(pc.y.v, testdatadict["ycoords"].v)
    assert np.allclose(pc.dirvec, testdatadict["dirvec"])


def test_export(testdatadict):
    """Test exporting the data to a dictionary."""
    pc = create_pointcloud(**testdatadict)
    dataout = pc.export()
    for key in testdatadict.keys():
        assert np.all(dataout[key] - testdatadict[key] == 0)


def test_init_hdf5(testdatadict):
    """Test initialization from HDF5."""
    def read_datasets_mock(qid, variant, key):
        match key:
            case "xcoords":
                return testdatadict["xcoords"]
            case "ycoords":
                return testdatadict["ycoords"]
            case "volume":
                return testdatadict["volume"]
            case "directionvector":
                return testdatadict["dirvec"]

    hdf5manager = MagicMock()
    hdf5manager.read_datasets.side_effect=read_datasets_mock
    pc = Pointcloud.init_hdf5()
    pc._treemanager = MagicMock()
    pc._treemanager.hdf5manager = hdf5manager

    assert all( np.isclose(pc.x, testdatadict["xcoords"], atol=0*unyt.m) )
    assert all( np.isclose(pc.y, testdatadict["ycoords"], atol=0*unyt.m) )
    assert all( np.isclose(pc.dirvec, testdatadict["dirvec"], atol=0) )
    assert np.isclose(pc.volume, testdatadict["volume"], atol=0*unyt.m**3)
    assert pc.npoint == testdatadict["xcoords"].size


def test_export_hdf5(testdatadict, hdf5manager):
    """Test exporting the data to a HDF5 file."""
    pc = create_pointcloud(**testdatadict)
    pc._treemanager = MagicMock()
    pc._treemanager.hdf5manager = hdf5manager
    pc._export_hdf5()
    pc_file = Pointcloud.init_hdf5()
    pc_file._treemanager = MagicMock()
    pc_file._treemanager.hdf5manager = hdf5manager

    assert all( np.isclose(pc.x, pc_file.x, atol=0*unyt.m) )
