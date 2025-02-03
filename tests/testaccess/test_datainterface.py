"""Test accessing data on disk or in C struct with unified interface."""
import os
import ctypes
from typing import Optional
from collections import namedtuple
from unittest.mock import MagicMock
from dataclasses import dataclass

import unyt
import pytest
import numpy as np

from a5py.exceptions import AscotIOException

from a5py import utils
from a5py.data.access import variants, InputVariant, Format
from a5py.data.access.hdf5 import HDF5Manager

from .conftest import QID1, DATE, NOTE, INPUTVAR, FNEMPTY

class Pointcloud(InputVariant):
    """A prototype class for testing purposes.

    All actual data classes inherit both from MetaDataHolder and DataHolder.
    Therefore, this class is a prototype for the data classes. This class
    doesn't have any real functionality beyond testing.
    """

    @dataclass
    class Struct(ctypes.Structure):
        """Data stored in C."""
        _pack_ = 1
        _fields_ = [
            ("Npoint", ctypes.c_int32),
            ("volume", ctypes.c_double),
            ("weight", ctypes.c_double),
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
    def weight(self):
        """Weight of the cloud."""
        if self._format == Format.HDF5:
            return self._read_hdf5("weight")
        return self._struct_.weight * unyt.kg

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
            "weight":self.weight,
        }
        return data


#pylint: disable=protected-access
def create_pointcloud(
        volume: utils.ArrayLike | None = None,
        dirvec: utils.ArrayLike | None = None,
        xcoords: utils.ArrayLike | None = None,
        ycoords: utils.ArrayLike | None = None,
        weight: Optional[utils.ArrayLike] = None,
        ) -> Pointcloud:
    """Creates a point cloud instance.

    This mimics the factory methods creating the actual data objects. Inputs
    should be checked in these methods.
    """
    parameters = variants.parse_parameters(
        volume, dirvec, xcoords, ycoords, weight,
        )
    ncoord = 1 if parameters["xcoords"] is None else parameters["xcoords"].size
    variants.validate_required_parameters(
            parameters,
            names=["volume", "dirvec", "xcoords", "ycoords"],
            units=["m**3", "1", "m", "m"],
            shape=[(), (3,), (ncoord,), (ncoord,)],
            dtype="f8",
            default=[
                0, np.array([1, 0, 0]), np.ones((1,)), np.ones((1,)),
                ],
        )
    variants.validate_optional_parameters(
        parameters, ["weight"], ["kg"], [()], "f8",
        [(parameters["volume"] * unyt.kg / unyt.m**3).v],
    )
    obj = Pointcloud()

    obj._struct_.Npoint = ncoord
    obj._struct_.volume = parameters["volume"]
    obj._struct_.weight = parameters["weight"]
    obj._struct_.dirvec[:] = parameters["dirvec"]

    offload_array_size = parameters["xcoords"].size + parameters["ycoords"].size
    obj._allocate_offload_arrays(offload_array_size, intsize=0)
    obj._append_offload_array(parameters["xcoords"], parameters["ycoords"])
    return obj


@pytest.fixture(name="testdata")
def fixture_testdata():
    """Test data for point cloud initialization."""
    Data = namedtuple(
        "data", ["xcoords", "ycoords", "dirvec", "volume", "weight"]
        )
    return Data(
        np.array([1.0, 2.0, 3.0, 4.0]),
        np.array([-1.0, -2.0, -3.0, -4.0]),
        np.array([1/np.sqrt(2), -1/np.sqrt(2), 0]),
        -5.0,
        2.0,
        )


@pytest.fixture(name="testdatadict")
def fixture_testdatadict(testdata):
    """Test data for creating point clouds in dictionary format + with units."""
    return {
        "volume":testdata.volume*unyt.m**3,
        "dirvec":testdata.dirvec,
        "xcoords":testdata.xcoords*unyt.m,
        "ycoords":testdata.ycoords*unyt.m,
        "weight":testdata.weight*unyt.kg,
        }


@pytest.fixture(name="pointcloud")
def fixture_pointcloud(testdata):
    """Create point cloud from test data."""
    pc = Pointcloud()

    pc._struct_.Npoint = testdata.xcoords.size
    pc._struct_.volume = testdata.volume
    pc._struct_.weight = testdata.weight
    pc._struct_.dirvec[:] = testdata.dirvec

    offload_array_size = testdata.xcoords.size + testdata.ycoords.size
    pc._allocate_offload_arrays(offload_array_size, intsize=0)
    pc._append_offload_array(testdata.xcoords, testdata.ycoords)
    return pc


@pytest.fixture(name="hdf5manager")
def fixture_hdf5manager():
    """Create HDF5 manager."""
    yield HDF5Manager(FNEMPTY, file_exists=False)
    os.unlink(FNEMPTY)


def test_init_and_access(pointcloud, testdata):
    """Test that initialized properties yield correct data."""
    pc = pointcloud
    assert pc.npoint == 4
    assert pc.volume.v == testdata.volume
    assert np.allclose(pc.x.v, testdata.xcoords)
    assert np.allclose(pc.y.v, testdata.ycoords)
    assert np.allclose(pc.dirvec, testdata.dirvec)


def test_struct_data_is_immutable(pointcloud):
    """Test that properties are immutable and return values that are copies."""
    value_to_ensure_data_is_copied = 10.0
    pc = pointcloud

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
    """Test the create method."""
    pc = create_pointcloud(**testdatadict)
    assert pc.npoint == 4
    assert np.allclose(pc.volume.v, testdatadict["volume"].v)
    assert np.allclose(pc.weight.v, testdatadict["weight"].v)
    assert np.allclose(pc.x.v, testdatadict["xcoords"].v)
    assert np.allclose(pc.y.v, testdatadict["ycoords"].v)
    assert np.allclose(pc.dirvec, testdatadict["dirvec"])


def test_create_different_but_acceptable_units(testdatadict):
    """Test the create method when using different units as expected."""
    volume = -2 * unyt.cm**3
    testdatadict["volume"] = volume
    pc = create_pointcloud(**testdatadict)
    assert np.allclose(pc.volume.v, testdatadict["volume"].to("m**3").v)


def test_create_not_acceptable_units(testdatadict):
    """Test the create method when using unaccetable units."""
    volume = -2 * unyt.m**2
    testdatadict["volume"] = volume
    with pytest.raises(ValueError):
        _ = create_pointcloud(**testdatadict)


@pytest.mark.parametrize(
        "coords", [[1,2,3], 1, np.array([[1,2,3]]), np.array([[1,2,3]]).T],
        ids=["list", "single element", "2D array", "2D array transposed"],
        )
def test_create_different_but_acceptable_shape(testdatadict, coords):
    """Test the create method when using different shape as expected."""
    testdatadict["xcoords"] = coords
    testdatadict["ycoords"] = coords
    pc = create_pointcloud(**testdatadict)
    assert np.allclose(pc.x.v, np.asarray(coords).ravel())


def test_create_inacceptable_or_inconsistent_shape(testdatadict):
    """Test the create method when parameters are inconsistent or incorrect."""
    testdatadict["xcoords"] = [1, 2, 3]
    testdatadict["ycoords"] = [1, 2]
    with pytest.raises(ValueError):
        _ = create_pointcloud(**testdatadict)

    testdatadict["ycoords"] = [1, 2, 3]
    testdatadict["volume"] = [1,2]
    with pytest.raises(ValueError):
        _ = create_pointcloud(**testdatadict)


def test_create_optional_missing(testdatadict):
    """Test that the optional parameter is properly set when not provided."""
    del testdatadict["weight"]
    pc = create_pointcloud(**testdatadict)
    assert np.allclose(pc.weight.v, testdatadict["volume"].v)


def test_create_required_missing(testdatadict):
    """Test that missing required parameter raises exception."""
    del testdatadict["volume"]
    with pytest.raises(ValueError):
        _ = create_pointcloud(**testdatadict)


def test_create_dummy():
    """Test the create method with dummy arguments."""
    pc = create_pointcloud()
    assert pc.npoint == 1
    assert np.allclose(pc.volume.v, 0.0)
    assert np.allclose(pc.weight.v, 0.0)
    assert np.allclose(pc.x.v, np.ones((1,)))
    assert np.allclose(pc.y.v, np.ones((1,)))
    assert np.allclose(pc.dirvec, np.array([1, 0, 0]))


def test_export(testdatadict):
    """Test exporting the data to a dictionary."""
    pc = create_pointcloud(**testdatadict)
    dataout = pc.export()
    for key in testdatadict.keys():
        assert np.all(dataout[key] - testdatadict[key] == 0)


def test_init_hdf5(testdatadict):
    """Test initialization from HDF5."""
    def read_datasets_mock(qid, variant, key):
        _, _ = qid, variant
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
    pc = Pointcloud()
    pc._format = Format.HDF5
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
    pc_file = Pointcloud()
    pc_file._format = Format.HDF5
    pc_file._treemanager = MagicMock()
    pc_file._treemanager.hdf5manager = hdf5manager

    assert all( np.isclose(pc.x, pc_file.x, atol=0*unyt.m) )
