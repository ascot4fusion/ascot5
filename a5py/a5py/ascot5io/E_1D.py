"""
Radial electric field IO.
"""
import h5py
import numpy as np
import random
import datetime
    
from . ascot5group import creategroup

def write_hdf5(fn, Nrho, r_eff, rhomin, rhomax, rho, dVdrho):
    """
    Write radial electric field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Nrho : int
        Number of rho slots in data.
    r_eff : real
        Number of rho slots in data.
    rhomin : real
        Minimum rho value.
    rhomax : real
        Maximum rho value.
    rho : real Nrho x 1 numpy array
        rho grid.
    dVdrho : real Nrho x 1 numpy array
        Gradient of electric potential in rho grid.
    """

    mastergroup = "efield"
    subgroup    = "E_1DS"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # TODO check that input is valid

    # Actual data.
    f[path].attrs['n_rho']   = Nrho
    f[path].attrs['r_eff']   = r_eff
    f[path].attrs['rho_min'] = rhomin
    f[path].attrs['rho_max'] = rhomax
    f.create_dataset(path + '/rho', data=rho, dtype='f8')
    f.create_dataset(path + '/dV_drho', data=dVdrho, dtype='f8')

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
    out["Nrho"]    = f[path].attrs['n_rho']
    out["r_eff"]   = f[path].attrs['r_eff']
    out["rho_min"] = f[path].attrs['rho_min']
    out["rho_max"] = f[path].attrs['rho_max']
    out["rho"]     = f[path + '/rho'][:]
    out["dVdrho"]  = f[path + '/dV_drho'][:]

    f.close()

    return out
