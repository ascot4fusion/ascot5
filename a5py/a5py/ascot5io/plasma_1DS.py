"""Plasma 1DS IO.
"""
import h5py
import numpy as np
import random
import datetime

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup


class plasma_1DS(DataGroup):
    """Object representing P_1DS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    @staticmethod
    def write_hdf5(fn, nrho, nion, anum, znum, mass, charge, rho,
               edensity, etemperature, idensity, itemperature, desc=None):
        """Write 1DS plasma input in HDF5 file.

        Parameters
        ----------
        fn : str
            Path to hdf5 file.
        nrho : int
            Number of rho grid points.
        nion : int
            Number of ion species.
        anum : array_like (nion,1)
            Ion species atomic mass number
        znum : array_like (nion,1)
            Ion species charge number.
        mass : array_like (nion,1)
            Ion species mass [amu].
        charge : array_like (nion,1)
            Ion species charge [e].
        rho : array_like (nrho,1)
            rho grid, doesn't have to be uniform.
        edensity : array_like (nrho,1)
            Electron density [m^-3].
        etemperature : array_like (nrho,1)
            Electron temperature [eV].
        idensity : array_like (nrho,nion)
            Ion density [m^-3].
        itemperature : array_like (nrho,1)
            Ion temperature [ev].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if etemperature.size != nrho:
            raise ValueError("Invalid size for electron temperature.")
        if itemperature.size != nrho:
            raise ValueError("Invalid size for ion temperature.")
        if edensity.size != nrho:
            raise ValueError("Invalid size for electron density.")
        if idensity.shape != (nrho,nion):
            raise ValueError("Invalid size for ion density.")

        idensity = np.transpose(idensity)

        parent = "plasma"
        group  = "plasma_1DS"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nion',   (1,1),    data=nion,   dtype='i4')
            g.create_dataset('nrho',   (1,1),    data=nrho,   dtype='i4')
            g.create_dataset('znum',   (nion,1), data=znum,   dtype='i4')
            g.create_dataset('anum',   (nion,1), data=anum,   dtype='i4')
            g.create_dataset('charge', (nion,1), data=charge, dtype='i4')
            g.create_dataset('mass',   (nion,1), data=mass,   dtype='f8')
            g.create_dataset('rho',    (nrho,1), data=rho,    dtype='f8')

            g.create_dataset('etemperature', (nrho,1),    data=etemperature,
                             dtype='f8')
            g.create_dataset('edensity',     (nrho,1),    data=edensity,
                             dtype='f8')
            g.create_dataset('itemperature', (nrho,1),    data=itemperature,
                             dtype='f8')
            g.create_dataset('idensity',     (nion,nrho), data=idensity,
                             dtype='f8')

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read P_1DS input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "plasma/plasma_1DS_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["idensity"] = np.transpose(out["idensity"])
        return out
