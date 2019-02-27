"""
Marker IO.

File: mrk_prt.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, n, ids, mass, charge,
               r, phi, z, vR, vphi, vz,
               anum, znum, weight, time, desc=None):
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
    group  = "prt"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("n",      (1,1), data=n,      dtype='i8').attrs['unit'] = '1';
        g.create_dataset("r",      (n,1), data=r,      dtype='f8').attrs['unit'] = 'm';
        g.create_dataset("phi",    (n,1), data=phi,    dtype='f8').attrs['unit'] = 'deg';
        g.create_dataset("z",      (n,1), data=z,      dtype='f8').attrs['unit'] = 'm';
        g.create_dataset("vr",     (n,1), data=vR,     dtype='f8').attrs['unit'] = 'm/s';
        g.create_dataset("vphi",   (n,1), data=vphi,   dtype='f8').attrs['unit'] = 'm/s';
        g.create_dataset("vz",     (n,1), data=vz,     dtype='f8').attrs['unit'] = 'm/s';
        g.create_dataset("mass",   (n,1), data=mass,   dtype='f8').attrs['unit'] = 'amu';
        g.create_dataset("charge", (n,1), data=charge, dtype='i4').attrs['unit'] = 'e';
        g.create_dataset("anum",   (n,1), data=anum,   dtype='i4').attrs['unit'] = '1';
        g.create_dataset("znum",   (n,1), data=znum,   dtype='i4').attrs['unit'] = '1';
        g.create_dataset("weight", (n,1), data=weight, dtype='f8').attrs['unit'] = 'markers/s';
        g.create_dataset("time",   (n,1), data=time,   dtype='f8').attrs['unit'] = 's';
        g.create_dataset("id",     (n,1), data=ids,    dtype='i8').attrs['unit'] = '1';


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

class mrk_prt(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
