"""
Radial electric field IO.

File: E_1D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, n_rho, rho_min, rho_max, dV_drho, r_eff, desc=None):
    """
    Write radial electric field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    n_rho : int
        Number of rho slots in data.
    r_eff : real
        Number of rho slots in data.
    rho_min : real
        Minimum rho value.
    rho_max : real
        Maximum rho value.
    rho : real n_rho x 1 numpy array
        rho grid.
    dV_drho : real n_rho x 1 numpy array
        Gradient of electric potential in rho grid.
    """

    parent = "efield"
    group  = "E_1DS"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc)

        g.create_dataset('n_rho',   data=n_rho,   dtype='i8')
        g.create_dataset('rho_min', data=rho_min, dtype='f8')
        g.create_dataset('rho_max', data=rho_max, dtype='f8')
        g.create_dataset('dV_drho', data=dV_drho, dtype='f8')
        g.create_dataset('r_eff',   data=r_eff,   dtype='f8')


def read_hdf5(fn, qid):
    """
    Read radial electric field input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the efield to be read.

    Returns
    -------

    Dictionary containing electric field data.
    """

    path = "efield" + "/E_1D-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["n_rho"]    = f[path + 'n_rho']
        out["rho_min"] = f[path + 'rho_min'][:]
        out["rho_max"] = f[path + 'rho_max'][:]
        out["dV_drho"]  = f[path + '/dV_drho'][:]
        out["r_eff"]   = f[path + 'r_eff'][:]

    return out

class E_1D(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
