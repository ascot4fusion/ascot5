"""
Spline radial electric field IO.

File: E_1DS.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, nrho, rhomin, rhomax, dvdrho, reff, desc=None):
    """
    Write radial electric field input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        nrho : int <br>
            Number of rho slots in data.
        rhomin : float <br>
            Minimum rho value.
        rhomax : float <br>
            Maximum rho value.
        dvdrho : array_like (nrho,1) <br>
            Gradient of electric potential in rho grid.
        reff : float <br>
            Number of rho slots in data.
        desc : str, optional <br>
            Input description.

    Returns:
        QID of the new input that was written.
    """

    parent = "efield"
    group  = "E_1DS"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc)

        g.create_dataset('nrho',   (1,),     data=nrho,   dtype='i8')
        g.create_dataset('rhomin', (1,),     data=rhomin, dtype='f8')
        g.create_dataset('rhomax', (1,),     data=rhomax, dtype='f8')
        g.create_dataset('dvdrho', (nrho,1), data=dvdrho, dtype='f8')
        g.create_dataset('reff',   (1,),     data=reff,   dtype='f8')


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


class E_1DS(AscotData):
    """
    Object representing E_1DS data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
