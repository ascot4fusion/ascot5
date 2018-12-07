"""
Marker IO.

File: mrk_fl.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.base import AscotInput

def write_hdf5(fn, N, ids, r, phi, z, pitch, weight, time, desc=None):
    """
    Write magnetic field line marker input in hdf5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    N : int
        Number of markers
    ids : int N x 1 numpy array
        unique identifier for each marker (positive integer)
    r : real N x 1 numpy array
        magnetic field line R coordinate
    phi : real N x 1 numpy array
        magnetic field line phi coordinate [deg]
    z : real N x 1 numpy array
        magnetic field line z coordinate
    pitch : real N x 1 numpy array
        magnetic field line pitch whose sign defines
        the direction field line is traced
    weight : real N x 1 numpy array
        magnetic field line weight (markers/s)
    time : real N x 1 numpy array
        magnetic field line initial time

    """
    parent = "marker"
    group  = "field_line"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group, desc=desc)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/n",      data=np.array([N]), dtype='i8').attrs['unit'] = '1';
        f.create_dataset(path + "/r",      data=r,             dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/phi",    data=phi,           dtype='f8').attrs['unit'] = 'deg';
        f.create_dataset(path + "/z",      data=z,             dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/pitch",  data=pitch,         dtype='f8').attrs['unit'] = '1';
        f.create_dataset(path + "/weight", data=weight,        dtype='f8').attrs['unit'] = 'markers/s';
        f.create_dataset(path + "/time",   data=time,          dtype='f8').attrs['unit'] = 's';
        f.create_dataset(path + "/id",     data=ids,           dtype='i8').attrs['unit'] = '1';

def read_hdf5(fn, qid):
    """
    Read field-line input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the field-line data to be read.

    Returns
    -------

    Dictionary containing field-line data.
    """

    out = {}
    with h5py.File(fn, "r") as f:
        path = "marker/field_line-"+qid

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for field in f[path]:
            out[field] = f[path][field][:]

    return out

class mrk_fl(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
