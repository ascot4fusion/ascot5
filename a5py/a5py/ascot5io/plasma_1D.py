"""
Plasma 1D IO.

File: plasma_1D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, nrho, nion, anum, znum, mass, charge, rho,
               edensity, etemperature, idensity, itemperature, desc=None):
    """
    Write 1D plasma input in HDF5 file.

    Args:
        fn : str <br>
            Path to hdf5 file.
        nrho : int <br>
            Number of rho grid points.
        nion : int <br>
            Number of ion species.
        anum : array_like (nion,1) <br>
            Ion species atomic mass number
        znum : array_like (nion,1) <br>
            Ion species charge number.
        mass : array_like (nion,1) <br>
            Ion species mass [amu].
        charge : array_like (nion,1) <br>
            Ion species charge [e].
        rho : array_like (nrho,1) <br>
            rho grid, doesn't have to be uniform.
        edensity : array_like (nrho,1) <br>
            Electron density [m^-3].
        etemperature : array_like (nrho,1) <br>
            Electron temperature [eV].
        idensity : array_like (nion,nrho) <br>
            Ion density [m^-3].
        itemperature : array_like (nrho,1) <br>
            Ion temperature [ev].
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "plasma"
    group  = "plasma_1D"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset('nion',   (1,1),    data=nion, dtype='i4')
        g.create_dataset('nrho',   (1,1),    data=nrho, dtype='i4')
        g.create_dataset('znum',   (Nion,1), data=znum, dtype='i4')
        g.create_dataset('anum',   (Nion,1), data=anum, dtype='i4')
        g.create_dataset('charge', (Nion,1), data=znum, dtype='i4')
        g.create_dataset('mass',   (Nion,1), data=anum, dtype='f8')
        g.create_dataset('rho',    (Nrho,1), data=rho,  dtype='f8')

        g.create_dataset('etemperature',   (nrho,1),    data=etemperature,
                         dtype='f8')
        g.create_dataset('edensity',       (nrho,1),    data=edensity,
                         dtype='f8')
        g.create_dataset('iontemperature', (nrho,1),    data=itemperature,
                         dtype='f8')
        g.create_dataset('iondensity',     (nion,nrho), data=idensity,
                         dtype='f8')

    return g.name


def read_hdf5(fn, qid):
    """
    Read P_1D input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "plasma/plasma_1D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out

class plasma_1D(AscotData):
    """
    Object representing P_1D data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
