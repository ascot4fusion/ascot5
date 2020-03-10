"""
3D wall IO.

File: wall_3D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

import a5py.wall.plot as plot

def write_hdf5(fn, nelements, x1x2x3, y1y2y3, z1z2z3, desc=None):
    """
    Write 3D wall input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        nelements : int <br>
            Number of wall triangles
        x1x2x3 : array_like (nelements,3) <br>
            Each triangle's vertices' x coordinates [m].
        y1y2y3 : array_like (nelements,3) <br>
            Each triangle's vertices' y coordinates [m].
        z1z2z3 : array_like (nelements,3) <br>
            Each triangle's vertices' z coordinates [m].
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """
    assert x1x2x3.shape == (nelements,3)
    assert y1y2y3.shape == (nelements,3)
    assert z1z2z3.shape == (nelements,3)

    parent = "wall"
    group  = "wall_3D"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset('x1x2x3',    (nelements,3), data=x1x2x3,    dtype='f8')
        g.create_dataset('y1y2y3',    (nelements,3), data=y1y2y3,    dtype='f8')
        g.create_dataset('z1z2z3',    (nelements,3), data=z1z2z3,    dtype='f8')
        g.create_dataset('nelements', (1,1),         data=nelements, dtype='i4')

    return gname


def read_hdf5(fn, qid):
    """
    Read 3D wall input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "wall/wall_3D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


class wall_3D(AscotData):
    """
    Object representing wall_3D data.
    """

    number_of_elements = None

    def read(self):
        W=read_hdf5(self._file, self.get_qid())
        self.number_of_elements =  W['nelements'][0][0]
        return W

    def getNumberOfElements(self):
        if self.number_of_elements is None:
            self.read()
        return self.number_of_elements


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)


    def area(self):
        w = self.read()

        ab_x = w["x1x2x3"][:,1] - w["x1x2x3"][:,0]
        ab_y = w["y1y2y3"][:,1] - w["y1y2y3"][:,0]
        ab_z = w["z1z2z3"][:,1] - w["z1z2z3"][:,0]

        ac_x = w["x1x2x3"][:,2] - w["x1x2x3"][:,0]
        ac_y = w["y1y2y3"][:,2] - w["y1y2y3"][:,0]
        ac_z = w["z1z2z3"][:,2] - w["z1z2z3"][:,0]

        A = 0.5 * np.sqrt(   (ab_y * ac_z - ab_z * ac_y)**2
                           + (ab_z * ac_x - ab_x * ac_z)**2
                           + (ab_x * ac_y - ab_y * ac_x)**2 )

        return A


    def plotRz(self, axes=None, phi=None):
        data = self.read()
        if phi is not None:
            plot.plot_intersection(data["x1x2x3"], data["y1y2y3"],
                                   data["z1z2z3"], phi, axes=axes)
        else:
            plot.plot_projection(data["x1x2x3"], data["y1y2y3"],
                                 data["z1z2z3"], axes=axes)
