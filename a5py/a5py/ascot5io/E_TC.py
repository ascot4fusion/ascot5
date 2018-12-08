"""
Trivial cartesian electric field IO.

The trivial cartesian electric field is a constant
electric field defined in cartesian coordinates.

File: E_TC.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotInput

def write_hdf5(fn, Exyz, desc=None):
    """
    Write trivial cartesian electric field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Exyz : real 3 x 1 numpy array
        Electric field value in cartesian coordinates
    """

    parent = "efield"
    group  = "E_TC"

    # Check that input is a valid array with three elements.
    if Exyz.shape != (3,1) and Exyz.shape != (1,3) and Exyz.shape != (3,):
        raise Exception('Exyz has invalid format.')

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("Exyz", (3,1), data = Exyz, dtype="f8")


def read_hdf5(fn, qid):
    """
    Read trivial cartesian electric field input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the efield to be read.

    Returns
    -------

    Dictionary containing electric field data.
    """

    path = "efield" + "/E_TC-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["Exyz"] = f[path]["Exyz"][:]

    return out

class E_TC(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
