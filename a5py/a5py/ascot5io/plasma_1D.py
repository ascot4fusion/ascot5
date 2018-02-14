"""
Plasma 1D IO.
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, Nrho, Nion, znum, anum, rho, ndens, ntemp, edens, etemp, idens, itemp):
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
    ndens : real Nrho x 1 numpy array
        neutral density (1/m^3) NOT IMPLEMNTED
    ntemp : real Nrho x 1 numpy array
        neutral temperature (eV) NOT IMPLEMENTED
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
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)
    
    # Neutrals are currently not implemented
    Nneutral = 1
    ndens = np.zeros((Nrho,Nneutral))
    ntemp = np.zeros((Nrho,1))
    
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
    f.create_dataset(path + '/n_neutrals', (1,1), dtype='i4', data=Nneutral)

    f.create_dataset(path + '/Z_num', (Nion,1), dtype='i4', data=znum)
    f.create_dataset(path + '/A_mass', (Nion,1), dtype='i4', data=anum)
    f.create_dataset(path + '/n_rho', (1,1), dtype='i4', data=Nrho)

    # 1D plasma properties
    f.create_dataset(path + '/rho', (Nrho,1), dtype='f8', data=rho)
    f.create_dataset(path + '/temp_0', (Nrho,1), dtype='f8', data=ntemp)
    f.create_dataset(path + '/dens_0', dtype='f8', data=ndens)
    f.create_dataset(path + '/temp_e', (Nrho,1), dtype='f8', data=etemp)
    f.create_dataset(path + '/dens_e', (Nrho,1), dtype='f8', data=edens)
    f.create_dataset(path + '/temp_i', (Nrho,1), dtype='f8', data=itemp)
    f.create_dataset(path + '/dens_i', dtype='f8', data=idens)
    
    f.close();


def read_hdf5(fn):
    """
    Read 1D plasma input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing plasma data.
    """

    group = "plasma"
    path = "plasma/P_1D"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

    # Actual data.
    out["Z_num"]    = f[group]["Z_num"][:]
    out["A_mass"]   = f[group]["A_mass"][:]
    out["Nion"]     = f[group].attrs['n_ions']
    out["Nneutral"] = f[group].attrs['n_neutrals']
    out["Nrho"]     = f[path].attrs['n_rho']

    out["rho"] = f[path]["rho"][:]
    out["ntemp"] = f[path]["temp_0"][:]
    out["ndens"] = f[path]["dens_0"][:]
    out["etemp"] = f[path]["temp_e"][:]
    out["edens"] = f[path]["dens_e"][:]
    out["itemp"] = f[path]["temp_i"][:]
    out["idens"] = np.transpose(np.reshape(f[path]["dens_i"][:], (out["Nion"], out["Nrho"]) ))

    f.close()

    return out
