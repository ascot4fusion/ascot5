"""
Marker IO.

File: mrk_prt.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.base import AscotInput

def write_hdf5(fn, N, ids, mass, charge,
               r, phi, z, vR, vphi, vz,
               weight, time, desc=None):
    """
    Write particle marker input in hdf5 file.

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
        particle R coordinate
    phi : real N x 1 numpy array
        particle phi coordinate [deg]
    z : real N x 1 numpy array
        particle z coordinate
    vR : real N x 1 numpy array
        particle velocity R-component
    vphi : real N x 1 numpy array
        particle velocity phi-component
    vz : real N x 1 numpy array
        particle velocity z-component
    weight : real N x 1 numpy array
        particle weight (markers/s)
    time : real N x 1 numpy array
        particle initial time

    """

    parent = "marker"
    group  = "particle"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group, desc=desc)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/n",      data=np.array([N]), dtype='i8').attrs['unit'] = '1';
        f.create_dataset(path + "/r",      data=r,             dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/phi",    data=phi,           dtype='f8').attrs['unit'] = 'deg';
        f.create_dataset(path + "/z",      data=z,             dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/v_r",    data=vR,            dtype='f8').attrs['unit'] = 'm/s';
        f.create_dataset(path + "/v_phi",  data=vphi,          dtype='f8').attrs['unit'] = 'm/s';
        f.create_dataset(path + "/v_z",    data=vz,            dtype='f8').attrs['unit'] = 'm/s';
        f.create_dataset(path + "/mass",   data=mass,          dtype='f8').attrs['unit'] = 'amu';
        f.create_dataset(path + "/charge", data=charge,        dtype='i4').attrs['unit'] = 'e';
        f.create_dataset(path + "/weight", data=weight,        dtype='f8').attrs['unit'] = 'markers/s';
        f.create_dataset(path + "/time",   data=time,          dtype='f8').attrs['unit'] = 's';
        f.create_dataset(path + "/id",     data=ids,           dtype='i8').attrs['unit'] = '1';

def read_hdf5(fn, qid):
    """
    Read particle input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the particle data to be read.

    Returns
    -------

    Dictionary containing particle data.
    """

    out = {}
    with h5py.File(fn, "r") as f:
        path = "marker/particle-"+qid

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for field in f[path]:
            out[field] = f[path][field][:]

    return out

class mrk_prt(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
