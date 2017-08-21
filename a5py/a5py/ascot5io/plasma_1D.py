"""
Plasma 1D IO.
"""
import h5py
import numpy as np
import random
import datetime

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
    
    # neutrals are currently not implemented
    ndens = np.zeros((Nrho,1))
    ntemp = np.zeros((Nrho,1))
    Nneutral = 1

    # check that input is valid
    if anum.size != Nion or znum.size != Nion:
        raise Exception('Number of ions in input not consistent')

    if rho.size != Nrho or edens.size != Nrho or etemp.size != Nrho or itemp.size != Nrho:
        raise Exception('Number of rho grid points in input not consistent')

    if Nrho*Nion != idens.size:
        raise Exception('Ion density data is not consisten with Nrho and Nion')

    if etemp[0] < 1 or etemp[0] > 1e5 or itemp[0] < 1 or itemp[0] >1e5:
        print "Warning: Check that temperature is given in eV"

    # convert ion density matrix in 1D array (which is how it is stored in hdf5)
    if idens.size != (Nion,Nrho):
        idens.flatten("C")
    else:
        idens.flatten("F")

    group = "plasma"
    type_ = "P_1D"
    path = "plasma/P_1D"
    
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

    # TODO Check that inputs are consistent.

    # Metadata.
    qid = random.getrandbits(64)
    f[path].attrs["qid"]  = np.int64_(qid)
    f[path].attrs["date"] = np.string_(datetime.datetime.now())


    f.create_dataset('plasma/Z_num', (Nion,1), dtype='i4', data=znum)
    f.create_dataset('plasma/A_mass', (Nion,1), dtype='i4', data=anum)
    f['plasma'].attrs['n_ions'] = Nion
    f['plasma'].attrs['n_neutrals'] = Nneutral

    # 1D plasma properties
    f.create_dataset(path + '/rho', (Nrho,1), dtype='f8', data=rho)
    f.create_dataset(path + '/temp_0', (Nrho,1), dtype='f8', data=ntemp)
    f.create_dataset(path + '/dens_0', (Nrho,1), dtype='f8', data=ndens)
    f.create_dataset(path + '/temp_e', (Nrho,1), dtype='f8', data=etemp)
    f.create_dataset(path + '/dens_e', (Nrho,1), dtype='f8', data=edens)
    f.create_dataset(path + '/temp_i', (Nrho,1), dtype='f8', data=itemp)
    f.create_dataset(path + '/dens_i', (Nrho*Nion,1), dtype='f8', data=idens)
    f[path].attrs['n_rho'] = Nrho

    f.close();


def read_hdf5(fn)
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
    out["Nion"]     = f['plasma'].attrs['n_ions']
    out["Nneutral"] = f['plasma'].attrs['n_neutrals']
    out["Nrho"]     = f[path].attrs['n_rho']

    out["rho"] = f[path]["rho"][:]
    out["ntemp"] = f[path]["temp_0"][:]
    out["ndens"] = f[path]["dens_0"][:]
    out["etemp"] = f[path]["temp_e"][:]
    out["edens"] = f[path]["dens_e"][:]
    out["itemp"] = f[path]["temp_i"][:]
    out["idens"] = f[path]["dens_i"][:]

    f.close()

    return out
