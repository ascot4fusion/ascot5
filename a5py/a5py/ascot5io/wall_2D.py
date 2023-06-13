"""2D wall IO.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

from . import wall_3D

import a5py.wall.plot as plot

class wall_2D(DataGroup):
    """
    Object representing wall_2D data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    def plotRz(self, axes=None, phi=0):
        w = self.read()
        plot.plot_segments(w["r"], w["z"], axes=axes)

    @staticmethod
    def write_hdf5(fn, nelements, r, z, desc=None):
        """Write 2D wall input in HDF5 file.

        First vertice shouldn't correspond to the last vertice, i.e., don't give
        (R,z) coordinates that make a closed loop.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        nelements : int
            Number of wall segments.
        r : array_like (n,1)
            R coordinates of wall segment vertices [m].
        z : array_like (n,1)
            z coordinates of wall segment vertices [m].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If the size of r or z does not match the number of elements.
        """
        if r.size != nelements or z.size != elements:
            raise ValueError("Size of r or z does not match the number of elements")

        parent = "wall"
        group  = "wall_2D"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("nelements", (1,1),         data=nelements, dtype='i4')
            g.create_dataset("r",         (nelements,1), data=r,         dtype='f8')
            g.create_dataset("z",         (nelements,1), data=z,         dtype='f8')

        return gname

    @staticmethod
    def write_hdf5_3D(fn, nelements, r, z, nphi, desc=None):

        x1x2x3 = np.ones((2*nelements*nphi, 3))
        y1y2y3 = np.ones((2*nelements*nphi, 3))
        r1r2r3 = np.ones((2*nelements*nphi, 3))
        p1p2p3 = np.ones((2*nelements*nphi, 3))
        z1z2z3 = np.ones((2*nelements*nphi, 3))

        def pol2cart(rho, phi):
            x = rho * np.cos(phi)
            y = rho * np.sin(phi)
            return(x, y)

        pv = np.linspace(0, 2*np.pi, nphi+1)
        for i in range(1,nphi+1):
            for j in range(n):
                r1r2r3[(i-1)*2*n + 2*(j-1),:] = [ r[j],    r[j],  r[j-1] ]
                p1p2p3[(i-1)*2*n + 2*(j-1),:] = [ pv[i-1], pv[i], pv[i] ]
                z1z2z3[(i-1)*2*n + 2*(j-1),:] = [ z[j],    z[j],  z[j-1] ]

                r1r2r3[(i-1)*2*n + 2*(j-1) + 1,:] = [ r[j],    r[j-1],  r[j-1] ]
                p1p2p3[(i-1)*2*n + 2*(j-1) + 1,:] = [ pv[i-1], pv[i-1], pv[i] ]
                z1z2z3[(i-1)*2*n + 2*(j-1) + 1,:] = [ z[j],    z[j-1],  z[j-1] ]

        x1x2x3,y1y2y3 = pol2cart(r1r2r3, p1p2p3)

        return {"nelements" : 2*nelements*nphi, "x1x2x3" : x1x2x3,
                "y1y2y3" : y1y2y3, "z1z2z3" : z1z2z3}

    @staticmethod
    def write_hdf5_dummy(fn, desc="Dummy"):
        r = np.array([0.01, 100, 100, 0.01])
        z = np.array([-100, -100, 100, 100])
        write_hdf5(fn=fn, nelements=4, r=r, z=z, desc=desc)
