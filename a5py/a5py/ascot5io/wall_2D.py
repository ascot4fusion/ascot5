"""
2D wall IO.

File: wall_2D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

from . import wall_3D

import a5py.wall.plot as plot

def write_hdf5(fn, nelements, r, z, desc=None):
    """
    Write 2D wall input in HDF5 file.

    First vertice shouldn't correspond to the last vertice, i.e., don't give
    (R,z) coordinates that make a closed loop.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        nelements : int <br>
            Number of wall segments.
        r : array_like (n,1) <br>
            R coordinates of wall segment vertices [m].
        z : array_like (n,1) <br>
            z coordinates of wall segment vertices [m].
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """
    assert r.size == nelements
    assert z.size == nelements

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


def read_hdf5(fn, qid):
    """
    Read 2D wall input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "wall/wall_2D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


def write_hdf5_3D(fn, n, R, z, nphi, desc=None):

    x1x2x3 = np.ones((2*n*nphi, 3))
    y1y2y3 = np.ones((2*n*nphi, 3))
    r1r2r3 = np.ones((2*n*nphi, 3))
    p1p2p3 = np.ones((2*n*nphi, 3))
    z1z2z3 = np.ones((2*n*nphi, 3))

    def pol2cart(rho, phi):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)

    pv = np.linspace(0, 2*np.pi, nphi+1)
    for i in range(1,nphi+1):
        for j in range(n):
            r1r2r3[(i-1)*2*n + 2*(j-1),:] = [ R[j],    R[j],  R[j-1] ]
            p1p2p3[(i-1)*2*n + 2*(j-1),:] = [ pv[i-1], pv[i], pv[i] ]
            z1z2z3[(i-1)*2*n + 2*(j-1),:] = [ z[j],    z[j],  z[j-1] ]

            r1r2r3[(i-1)*2*n + 2*(j-1) + 1,:] = [ R[j],    R[j-1],  R[j-1] ]
            p1p2p3[(i-1)*2*n + 2*(j-1) + 1,:] = [ pv[i-1], pv[i-1], pv[i] ]
            z1z2z3[(i-1)*2*n + 2*(j-1) + 1,:] = [ z[j],    z[j-1],  z[j-1] ]

    x1x2x3,y1y2y3 = pol2cart(r1r2r3, p1p2p3)
    wall_3D.write_hdf5(fn, 2*n*nphi, x1x2x3, y1y2y3, z1z2z3, desc=desc)


class wall_2D(AscotData):
    """
    Object representing wall_2D data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())

    def plotRz(self, axes=None, phi=0):
        w = self.read()
        plot.plot_segments(w["r"], w["z"], axes=axes)
