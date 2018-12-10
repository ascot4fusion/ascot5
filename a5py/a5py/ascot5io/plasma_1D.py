"""
Plasma 1D IO.

File: plasma_1D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, Nrho, Nion, znum, anum, rho, edens, etemp, idens, itemp,
               desc=None):
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

    parent = "plasma"
    group  = "plasma_1D"


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


    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset('n_ions', (1,1),    data=Nion,  dtype='i4')
        g.create_dataset('Z_num',  (Nion,1), data=znum,  dtype='i4')
        g.create_dataset('A_mass', (Nion,1), data=anum,  dtype='i4')
        g.create_dataset('n_rho',  (1,1),    data=Nrho,  dtype='i4')
        g.create_dataset('rho',    (Nrho,1), data=rho,   dtype='f8')
        g.create_dataset('temp_e', (Nrho,1), data=etemp, dtype='f8')
        g.create_dataset('dens_e', (Nrho,1), data=edens, dtype='f8')
        g.create_dataset('temp_i', (Nrho,1), data=itemp, dtype='f8')
        g.create_dataset('dens_i',           data=idens, dtype='f8')


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

class plasma_1D(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
