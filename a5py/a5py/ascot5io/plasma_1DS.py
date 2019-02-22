"""
Plasma 1DS IO.

File: plasma_1DS.py
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, Nrho, Nion, znum, anum, rhomin, rhomax,
               edens, etemp, idens, itemp, desc=None):
    """
    Write 1DS plasma input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        Nrho : int <br>
            Number of rho grid points.
        Nion : int <br>
            Number of ions.
        znum : int, array_like (Nion,) <br>
            Ion charge numbers.
        anum : int, array_like (Nion,) <br>
            Ion mass numbers.
        rhomin : real <br>
            Minimum rho grid edge value.
        rhomax : real <br>
            Maximum rho grid edge value.
        edens : real array_like (Nrho,)
            Electron density [1/m^3].
        etemp : real array_like (Nrho,)
            Electron temperature [eV].
        idens : real array_like (Nrho, Nion)
            Ion density [1/m^3].
        itemp : real array_like (Nrho,)
            Ion temperature [eV]
    """

    parent = "plasma"
    group  = "plasma_1DS"

    #idens = np.transpose(idens)

    if etemp[0] < 1 or etemp[0] > 1e5 or itemp[0] < 1 or itemp[0] >1e5:
        print("Warning: Check that temperature is given in eV")


    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset('n_ions',  (1,1),       data=Nion,   dtype='i4')
        g.create_dataset('Z_num',   (Nion,1),    data=znum,   dtype='i4')
        g.create_dataset('A_mass',  (Nion,1),    data=anum,   dtype='i4')
        g.create_dataset('n_rho',   (1,1),       data=Nrho,   dtype='i4')
        g.create_dataset('rho_min', (1,1),       data=rhomin, dtype='f8')
        g.create_dataset('rho_max', (1,1),       data=rhomax, dtype='f8')
        g.create_dataset('temp_e',  (Nrho,1),    data=etemp,  dtype='f8')
        g.create_dataset('dens_e',  (Nrho,1),    data=edens,  dtype='f8')
        g.create_dataset('temp_i',  (Nrho,1),    data=itemp,  dtype='f8')
        g.create_dataset('dens_i',  (Nrho,Nion), data=idens,  dtype='f8')


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


class plasma_1DS(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
