"""
Marker IO.

File: mrk_gc.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.base import AscotInput

def write_hdf5(fn, N, ids, mass, charge,
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

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group, desc=desc)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/n",      data=np.array([N]), dtype='i8').attrs['unit'] = '1';
        f.create_dataset(path + "/r",      data=r,             dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/phi",    data=phi,           dtype='f8').attrs['unit'] = 'deg';
        f.create_dataset(path + "/z",      data=z,             dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/energy", data=energy,        dtype='f8').attrs['unit'] = 'ev';
        f.create_dataset(path + "/pitch",  data=pitch,         dtype='f8').attrs['unit'] = '1';
        f.create_dataset(path + "/theta",  data=theta,         dtype='f8').attrs['unit'] = 'rad';
        f.create_dataset(path + "/mass",   data=mass,          dtype='f8').attrs['unit'] = 'amu';
        f.create_dataset(path + "/charge", data=charge,        dtype='i4').attrs['unit'] = 'e';
        f.create_dataset(path + "/weight", data=weight,        dtype='f8').attrs['unit'] = 'markers/s';
        f.create_dataset(path + "/time",   data=time,          dtype='f8').attrs['unit'] = 's';
        f.create_dataset(path + "/id",     data=ids,           dtype='i8').attrs['unit'] = '1';

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
