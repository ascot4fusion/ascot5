"""
Radial electric field IO.
"""
import h5py
import numpy as np
import random
import datetime
    
def write_hdf5(fn, Nrho, r_eff, rhomin, rhomax, rho, dVdrho):
    """
    Write radial electric field input in HDF5 file.

    TODO Not compatible with new HDF5 format.

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

    group = "efield"
    type_ = "E_1D"
    path = "efield/E_1D"

    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    if not group in f:
        o = f.create_group(group)
        o.attrs["type"] = np.string_(type_)
    else:
        o = f[group]
        del o.attrs["type"]
        o.attrs["type"] = np.string_(type_)

    # Remove group if one is already present.
    if path in f:
        del f[path]
    f.create_group(path)

    # TODO check that input is valid

    # Metadata.
    qid = random.getrandbits(64)
    f[path].attrs["qid"]  = np.int64_(qid)
    f[path].attrs["date"] = np.string_(datetime.datetime.now())

    # Actual data.
    f[path].attrs['n_rho']   = Nrho
    f[path].attrs['r_eff']   = r_eff
    f[path].attrs['rho_min'] = rhomin
    f[path].attrs['rho_max'] = rhomax
    f.create_dataset(path + '/rho', data=rho, dtype='f8')
    f.create_dataset(path + '/dV_drho', data=dVdrho, dtype='f8')

    f.close()


def read_hdf5(fn):
    """
    Read radial electric field input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing electric field data.
    """

    path = "efield/E_1D"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]
    
    # Actual data.
    out["Nrho"]    = f[path].attrs['n_rho']
    out["r_eff"]   = f[path].attrs['r_eff']
    out["rho_min"] = f[path].attrs['rho_min']
    out["rho_max"] = f[path].attrs['rho_max']
    out["rho"]     = f[path + '/rho'][:]
    out["dVdrho"]  = f[path + '/dV_drho'][:]

    f.close()

    return out
