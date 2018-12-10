"""
Orbits HDF5 IO module.

File: orbits.py
"""
import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid, name):
    """
    Read orbits.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    qid : str
        qid of the run these orbit data correspond to.

    Returns
    -------

    Dictionary storing the orbits that were read.
    """

    with h5py.File(fn,"r") as f:
        orbits = f["/results/run-"+qid+"/orbits/" + name]

        out = {}

        # Read data from file.
        for field in orbits:
            out[field]           = orbits[field][:]
            out[field + "_unit"] = orbits[field].attrs["unit"]

        # Find how many markers we have and their ids.
        out["N"]        = (np.unique(orbits["id"][:])).size
        out["uniqueId"] = np.unique(orbits["id"][:])

        # Sort fields by id (major) and time (minor), both ascending.
        if out["N"] > 0:
            ind = np.lexsort((out["time"], out["id"]))
            for field in orbits:
                if field[-4:] != "unit":
                    out[field] = out[field][ind]

    return out


def write_hdf5(fn, orbits, qid):
    """
    Write orbits.

    Unlike most other "write" functions, this one takes dictionary
    as an argument. The dictionary should have exactly the same format
    as given by the "read" function in this module. The reason for this
    is that this function is intended to be used only when combining
    different HDF5 files into one.

    TODO not compatible with new format

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    orbits : dictionary
        Orbit data to be written in dictionary format.
    qid : int
        Run id these orbits correspond to.
    """

    with h5py.File(fn, "a") as f:

        for orbgroup in orbits:
            path = "results/run-" + qid + "/orbits/" + orbgroup

            # Remove group if one is already present.
            if path in f:
                del f[path]
            f.create_group(path)

            # TODO Check that inputs are consistent.

            # Write data to file.
            for field in orbits[orbgroup]:
                if field[-4:] != "unit" and field != "N" and field != "uniqueId":
                    d = f.create_dataset(path + "/" + field, data=orbits[orbgroup][field])
                    d.attrs["unit"] = orbits[orbgroup][field + "_unit"]

class Orbits(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid(), self.get_type())
