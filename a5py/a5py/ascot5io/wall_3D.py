"""
3D wall IO.

File: wall_3D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, n, x1x2x3, y1y2y3, z1z2z3, flag, desc=None):
    """
    Write 3D wall input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    n : int
        Number of wall triangles
    x1x2x3, y1y2y3, z1z2z3 : real n x 3 numpy array
        Wall triangle vertices' xyz coordinates.
    flag : int n x 1 numpy array
        Indicates which part of the wall (e.g. divertor) triangle belongs to.
    """

    parent = "wall"
    group  = "wall_3D"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset('x1x2x3',    (n,3), data=x1x2x3,          dtype='f8')
        g.create_dataset('y1y2y3',    (n,3), data=y1y2y3,          dtype='f8')
        g.create_dataset('z1z2z3',    (n,3), data=z1z2z3,          dtype='f8')
        g.create_dataset('flag',      (n,1), data=flag,            dtype='i4')
        g.create_dataset('nelements', (1,),  data=n,               dtype='i4')
        g.create_dataset('xmin',      (1,),  data=np.amin(x1x2x3), dtype='f8')
        g.create_dataset('xmax',      (1,),  data=np.amax(x1x2x3), dtype='f8')
        g.create_dataset('ymin',      (1,),  data=np.amin(y1y2y3), dtype='f8')
        g.create_dataset('ymax',      (1,),  data=np.amax(y1y2y3), dtype='f8')
        g.create_dataset('zmin',      (1,),  data=np.amin(z1z2z3), dtype='f8')
        g.create_dataset('zmax',      (1,),  data=np.amax(z1z2z3), dtype='f8')


def read_hdf5(fn, qid):
    """
    Read 3D wall input from HDF5 file.

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

    path = "wall" + "/wall_3D-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["n"]      = f[path]["n"][:]
        out["x1x2x3"] = f[path]["x1x2x3"][:]
        out["y1y2y3"] = f[path]["y1y2y3"][:]
        out["z1z2z3"] = f[path]["z1z2z3"][:]
        out["flag"]   = f[path]["flag"][:]

    return out

class wall_3D(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
