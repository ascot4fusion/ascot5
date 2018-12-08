"""
Marker IO.

File: mrk_fl.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotInput

def write_hdf5(fn, n, ids, r, phi, z, pitch, weight, time, desc=None):
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

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("n",      (1,1), data=n,      dtype='i8').attrs['unit'] = '1';
        g.create_dataset("r",             data=r,      dtype='f8').attrs['unit'] = 'm';
        g.create_dataset("phi",           data=phi,    dtype='f8').attrs['unit'] = 'deg';
        g.create_dataset("z",             data=z,      dtype='f8').attrs['unit'] = 'm';
        g.create_dataset("pitch",         data=pitch,  dtype='f8').attrs['unit'] = '1';
        g.create_dataset("weight",        data=weight, dtype='f8').attrs['unit'] = 'markers/s';
        g.create_dataset("time",          data=time,   dtype='f8').attrs['unit'] = 's';
        g.create_dataset("id",            data=ids,    dtype='i8').attrs['unit'] = '1';

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
