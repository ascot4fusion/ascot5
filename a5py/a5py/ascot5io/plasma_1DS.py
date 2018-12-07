"""
Plasma 1DS IO.

File: plasma_1DS.py
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5file import add_group

def write_hdf5(fn, Nrho, Nion, znum, anum, rhomin, rhomax,
               ndens, ntemp, edens, etemp, idens, itemp):
    """
    Write 1DS plasma input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Nrho : int
        Number of rho grid points.
    Nion : int
        Number of ions.
    znum : int Nion x 1 numpy array
        Ion charge numbers.
    anum : int Nion x 1 numpy array
        Ion mass numbers.
    rhomin, rhomax : real
        Edges of the uniform rho-grid.
    ndens : real Nrho x 1 numpy array
        Neutral density (1/m^3) NOT IMPLEMNTED
    ntemp : real Nrho x 1 numpy array
        Neutral temperature (eV) NOT IMPLEMENTED
    edens : real Nrho x 1 numpy array
        Electron density (1/m^3).
    etemp : real Nrho x 1 numpy array
        Electron temperature (eV).
    idens : real Nrho x Nion numpy array
        Ion density (1/m^3).
    itemp : real Nrho x 1 numpy array
        Ion temperature (eV)
    """

    parent = "plasma"
    group  = "plasma_1DS"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group)

        # Neutrals are currently not implemented
        Nneutral = 1
        ndens = np.zeros((Nrho,Nneutral))
        ntemp = np.zeros((Nrho,1))

        # Check that input is valid
        if anum.size != Nion or znum.size != Nion:
            raise Exception('Number of ions in input not consistent')

        if (rho.size != Nrho or edens.size != Nrho
            or etemp.size != Nrho or itemp.size != Nrho):
            raise Exception('Number of rho grid points in input not consistent')

        if Nrho != idens.shape[0] or Nion != idens.shape[1]:
            idens = np.transpose(idens)
            if Nrho != idens.shape[0] or Nion != idens.shape[1]:
                raise Exception('Ion density data is not consistent'
                                'with Nrho and Nion')

        idens = np.transpose(idens)

        if etemp[0] < 1 or etemp[0] > 1e5 or itemp[0] < 1 or itemp[0] >1e5:
            print("Warning: Check that temperature is given in eV")


        # TODO Check that inputs are consistent.

        f.create_dataset(path + '/n_ions', (1,1), dtype='i4', data=Nion)
        f.create_dataset(path + '/n_neutrals', (1,1), dtype='i4', data=Nneutral)

        f.create_dataset(path + '/Z_num', (Nion,1), dtype='i4', data=znum)
        f.create_dataset(path + '/A_mass', (Nion,1), dtype='i4', data=anum)
        f.create_dataset(path + '/n_rho', (1,1), dtype='i4', data=Nrho)
        f.create_dataset(path + '/rho_min', (1,1), dtype='f8', data=rhomin)
        f.create_dataset(path + '/rho_max', (1,1), dtype='f8', data=rhomax)

        # 1DS plasma properties
        f.create_dataset(path + '/temp_0', (Nrho,1), dtype='f8', data=ntemp)
        f.create_dataset(path + '/dens_0', dtype='f8', data=ndens)
        f.create_dataset(path + '/temp_e', (Nrho,1), dtype='f8', data=etemp)
        f.create_dataset(path + '/dens_e', (Nrho,1), dtype='f8', data=edens)
        f.create_dataset(path + '/temp_i', (Nrho,1), dtype='f8', data=itemp)
        f.create_dataset(path + '/dens_i', dtype='f8', data=idens)


def read_hdf5(fn, qid):
    """
    Read 1DS plasma input from HDF5 file.

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

    path = "plasma" + "/plasma_1DS-" + qid

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
        out["Nneutral"] = f[path]['n_neutrals'][:]
        out["Nrho"]     = f[path]['n_rho'][:]
        out["rhomin"]     = f[path]['rhomin'][:]
        out["rhomax"]     = f[path]['rhomax'][:]

        out["ntemp"] = f[path]["temp_0"][:]
        out["ndens"] = f[path]["dens_0"][:]
        out["etemp"] = f[path]["temp_e"][:]
        out["edens"] = f[path]["dens_e"][:]
        out["itemp"] = f[path]["temp_i"][:]
        out["idens"] = f[path]["dens_i"][:]

    return out
