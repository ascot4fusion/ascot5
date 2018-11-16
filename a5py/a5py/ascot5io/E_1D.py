"""
Radial electric field IO.
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, n_rho, rho_min, rho_max, dV_drho, r_eff):
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

    mastergroup = "efield"
    subgroup    = "E_1DS"

    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # TODO check that input is valid

    # Actual data.
    f.create_dataset(path + '/n_rho', data=n_rho, dtype='i8')
    f.create_dataset(path + '/rho_min', data=rho_min, dtype='f8')
    f.create_dataset(path + '/rho_max', data=rho_max, dtype='f8')
    f.create_dataset(path + '/dV_drho', data=dV_drho, dtype='f8')
    f.create_dataset(path + '/r_eff', data=r_eff, dtype='f8')

    f.close()


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

    f = h5py.File(fn,"r")

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

    f.close()

    return out
