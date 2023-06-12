"""Spline radial electric field IO.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from . E import E

def write_hdf5(fn, nrho, rhomin, rhomax, dvdrho, reff, desc=None):
    """Write radial electric field input in HDF5 file.

    Parameters
    ----------
    fn : str
        Full path to the HDF5 file.
    nrho : int
        Number of rho slots in data.
    rhomin : float
        Minimum rho value.
    rhomax : float
        Maximum rho value.
    dvdr : array_like (nrho,1)
        Derivative of electric potential WRT minor radius [V/m].

        If reff = 1, this is essentially equal to dv/drho
    reff : float
        Effective minor radius of the plasma [m].
    desc : str, optional
        Input description.

    Returns
    -------
    name : str
        Name, i.e. "<type>_<qid>", of the new input that was written.
    """

    parent = "efield"
    group  = "E_1DS"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc)
        gname = g.name.split("/")[-1]

        g.create_dataset('nrho',   (1,),     data=nrho,   dtype='i8')
        g.create_dataset('rhomin', (1,),     data=rhomin, dtype='f8')
        g.create_dataset('rhomax', (1,),     data=rhomax, dtype='f8')
        g.create_dataset('dvdrho', (nrho,1), data=dvdrho, dtype='f8')
        g.create_dataset('reff',   (1,),     data=reff,   dtype='f8')

    return gname


def read_hdf5(fn, qid):
    """
    Read radial electric field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "efield/E_1DS_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


class E_1DS(E):
    """
    Object representing E_1DS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
