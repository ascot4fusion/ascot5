"""
Marker IO.

File: mrk_gc.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotInput

def write_hdf5(fn, n, ids, mass, charge,
               r, phi, z, energy, pitch, theta,
               weight, time, desc=None):
    """
    Write guiding center marker input in hdf5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    N : int
        Number of markers
    ids : int N x 1 numpy array
        unique identifier for each marker (positive integer)
    charge : int
        charge (e)
    mass : real
        mass (amu)
    r : real N x 1 numpy array
        guiding center R coordinate
    phi : real N x 1 numpy array
        guiding center phi coordinate [deg]
    z : real N x 1 numpy array
        guiding center z coordinate
    energy : real N x 1 numpy array
        guiding center energy (eV)
    pitch : real N x 1 numpy array
        guiding center pitch (v_para/v_tot)
    theta : real N x 1 numpy array
        guiding center gyroangle (rad)
    weight : real N x 1 numpy array
        guiding center weight (markers/s)
    time : real N x 1 numpy array
        guiding center initial time

    """
    parent = "marker"
    group  = "guiding_center"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("n",      (1,1), data=n,      dtype='i8').attrs['unit'] = '1';
        g.create_dataset("r",             data=r,      dtype='f8').attrs['unit'] = 'm';
        g.create_dataset("phi",           data=phi,    dtype='f8').attrs['unit'] = 'deg';
        g.create_dataset("z",             data=z,      dtype='f8').attrs['unit'] = 'm';
        g.create_dataset("energy",        data=energy, dtype='f8').attrs['unit'] = 'ev';
        g.create_dataset("pitch",         data=pitch,  dtype='f8').attrs['unit'] = '1';
        g.create_dataset("theta",         data=theta,  dtype='f8').attrs['unit'] = 'rad';
        g.create_dataset("mass",          data=mass,   dtype='f8').attrs['unit'] = 'amu';
        g.create_dataset("charge",        data=charge, dtype='i4').attrs['unit'] = 'e';
        g.create_dataset("weight",        data=weight, dtype='f8').attrs['unit'] = 'markers/s';
        g.create_dataset("time",          data=time,   dtype='f8').attrs['unit'] = 's';
        g.create_dataset("id",            data=ids,    dtype='i8').attrs['unit'] = '1';

def read_hdf5(fn, qid):
    """
    Read guiding-center input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the guiding-center data to be read.

    Returns
    -------

    Dictionary containing guiding-center data.
    """

    out = {}
    with h5py.File(fn, "r") as f:
        path = "marker/guiding_center-"+qid

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for field in f[path]:
            out[field] = f[path][field][:]

    return out

class mrk_gc(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
