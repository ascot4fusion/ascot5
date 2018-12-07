"""
2D wall IO.

File: wall_2D.py
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5file import add_group
from a5py.ascot5io.base import AscotInput

from . import wall_3D

def write_hdf5(fn, n, R, z, desc=None):
    """
    Write 2D wall input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    n : int
        Number of wall segments
    R, z : real n x 1 numpy array
        Wall segment vertices' R and z coordinates

    Notes
    -----

    First (R,z) point doesn't have to correspond to last point
    as the wall is closed automatically. TODO check this
    """

    parent = "wall"
    group  = "wall_2D"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group, desc=desc)

        # TODO Check that inputs are consistent.

        f.create_dataset(path + "/n", (1,1), data=n, dtype='i4')
        f.create_dataset(path + "/r", data=R, dtype='f8')
        f.create_dataset(path + "/z", data=z, dtype='f8')

        setdescription(f, mastergroup, desc)


def read_hdf5(fn, qid):
    """
    Read 2D wall input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the wall to be read.

    Returns
    -------

    Dictionary containing wall data.
    """

    path = "wall" + "/wall_2D-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["n"] = f[path]["n"][:]
        out["R"] = f[path]["r"][:]
        out["z"] = f[path]["z"][:]

    return out

def write_hdf5_3D(fn, n, R, z, nphi, desc=None):

    x1x2x3 = np.ones((2*n*nphi, 3))
    y1y2y3 = np.ones((2*n*nphi, 3))
    r1r2r3 = np.ones((2*n*nphi, 3))
    p1p2p3 = np.ones((2*n*nphi, 3))
    z1z2z3 = np.ones((2*n*nphi, 3))
    flag   = np.ones((2*n*nphi, 1))

    def pol2cart(rho, phi):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)

    pv = np.linspace(0, 2*np.pi, nphi+1)
    for i in range(1,nphi):
        for j in range(1,n-1):
            r1r2r3[i-1 + 2*(j-1),:] = [ R[j], R[j], R[j-1] ]
            p1p2p3[i-1 + 2*(j-1),:] = [ pv[i-1], pv[i], pv[i] ]
            z1z2z3[i-1 + 2*(j-1),:] = [ z[j], z[j], z[j-1] ]

            r1r2r3[i-1 + 2*(j-1) + 1,:] = [ R[j], R[j-1], R[j-1] ]
            p1p2p3[i-1 + 2*(j-1) + 1,:] = [ pv[i-1], pv[i-1], pv[i] ]
            z1z2z3[i-1 + 2*(j-1) + 1,:] = [ z[j], z[j-1], z[j-1] ]

    x1x2x3,y1y2y3 = pol2cart(p1p2p3, r1r2r3)
    wall_3D.write_hdf5(fn, 2*n*nphi, x1x2x3, y1y2y3, z1z2z3, flag, desc=desc)

class wall_2D(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
