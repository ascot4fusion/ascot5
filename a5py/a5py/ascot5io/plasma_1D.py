"""
Plasma 1D IO.
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, Nrho, Nion, znum, anum, rho, edens, etemp, idens, itemp):
    """
    Write 1D plasma input in HDF5 file.

    Parameters
    ----------

    fn : str
        path to hdf5 file
    Nrho : int
        Number of rho grid points
    Nion : int
        Number of ions
    znum : int Nion x 1 numpy array
        Ion charge number
    anum : int Nion x 1 numpy array
        Ion mass number
    rho : int Nrho x 1 numpy array
        rho grid array
    edens : real Nrho x 1 numpy array
        electron density (1/m^3)
    etemp : real Nrho x 1 numpy array
        electron temperature (eV)
    idens : real Nrho x Nion numpy array
        ion density (1/m^3)
    itemp : real Nrho x 1 numpy array
        ion temperature (eV)
    """

    mastergroup = "plasma"
    subgroup    = "plasma_1D"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = creategroup(f, mastergroup, subgroup)

        # Check that input is valid
        if anum.size != Nion or znum.size != Nion:
            raise Exception('Number of ions in input not consistent')

        if rho.size != Nrho or edens.size != Nrho or etemp.size != Nrho or itemp.size != Nrho:
            raise Exception('Number of rho grid points in input not consistent')

        if Nrho != idens.shape[0] or Nion != idens.shape[1]:
            idens = np.transpose(idens)
            if Nrho != idens.shape[0] or Nion != idens.shape[1]:
                raise Exception('Ion density data is not consisten with Nrho and Nion')

        idens = np.transpose(idens)

        if etemp[0] < 1 or etemp[0] > 1e5 or itemp[0] < 1 or itemp[0] >1e5:
            print("Warning: Check that temperature is given in eV")


        # TODO Check that inputs are consistent.

        f.create_dataset(path + '/n_ions', (1,1), dtype='i4', data=Nion)

        f.create_dataset(path + '/Z_num', (Nion,1), dtype='i4', data=znum)
        f.create_dataset(path + '/A_mass', (Nion,1), dtype='i4', data=anum)
        f.create_dataset(path + '/n_rho', (1,1), dtype='i4', data=Nrho)

        # 1D plasma properties
        f.create_dataset(path + '/rho', (Nrho,1), dtype='f8', data=rho)
        f.create_dataset(path + '/temp_e', (Nrho,1), dtype='f8', data=etemp)
        f.create_dataset(path + '/dens_e', (Nrho,1), dtype='f8', data=edens)
        f.create_dataset(path + '/temp_i', (Nrho,1), dtype='f8', data=itemp)
        f.create_dataset(path + '/dens_i', dtype='f8', data=idens)


def read_hdf5(fn, qid):
    """
    Read 1D plasma input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the plasma to be read.

    Returns
    -------

    Dictionary containing plasma data.
    """

    path = "plasma" + "/plasma_1D-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["Z_num"]    = f[path]["Z_num"][:]
        out["A_mass"]   = f[path]["A_mass"][:]
        out["Nion"]     = f[path]['n_ions'][:]
        out["Nrho"]     = f[path]['n_rho'][:]

        out["rho"] = f[path]["rho"][:]
        out["etemp"] = f[path]["temp_e"][:]
        out["edens"] = f[path]["dens_e"][:]
        out["itemp"] = f[path]["temp_i"][:]
        out["idens"] = f[path]["dens_i"][:]

    return out
