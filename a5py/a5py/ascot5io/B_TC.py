"""Trivial Cartesian magnetic field HDF5 IO.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

class B_TC(DataGroup):
    """
    Object representing B_TC data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    @staticmethod
    def write_hdf5(fn, bxyz, jacobian, rhoval, psival=None, axisr=1, axisz=0,
                   desc=None):
        """Write trivial cartesian magnetic field input to HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        bxyz : array_like (3,1)
            Magnetic field in cartesian coordinates at origo.
        jacobian : array_like (3,3)
            Magnetic field Jacobian, jacobian[i,j] = dB_i/dx_j
        rhoval: float
            Constant rho value.
        psival: float, optional
            Constant psi value. If None, same as rhoval.
        axisr: float, optional
            Magnetic axis R coordinate.
        axisz: real, optional
            Magnetic axis z coordinate.
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
        if bxyz.shape != (3,1): raise ValueError("Invalid shape for Bxyz.")

        parent = "bfield"
        group  = "B_TC"
        gname  = ""

        if psival is None:
            psival = rhoval

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("bxyz",     (3,1), data=bxyz,     dtype="f8")
            g.create_dataset("jacobian", (3,3), data=jacobian, dtype="f8")
            g.create_dataset("rhoval",   (1,),  data=rhoval,   dtype="f8")
            g.create_dataset("psival",   (1,),  data=psival,   dtype="f8")
            g.create_dataset("axisr",    (1,),  data=axisr,    dtype="f8")
            g.create_dataset("axisz",    (1,),  data=axisz,    dtype="f8")

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read Cartesian magnetic field input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "bfield/B_TC_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        return out

    @staticmethod
    def write_hdf5_dummy(fn, desc="Dummy"):
        """
        Write dummy data.
        """
        B   = np.array([1,0,3])
        jac = np.array([ [0,0,0.01], [0,0,0], [0,0,0] ])
        return write_hdf5(fn, B, jacobian=jac,
                          rhoval=0.5, psival=1.5, axisr=6, axisz=0, desc=desc)
