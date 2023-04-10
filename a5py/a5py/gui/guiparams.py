"""
Store and retrieve GUI related parameters to/from HDF5
"""

import numpy as np
import h5py

def store(hdf5, **kwargs):
    """
    Stores parameter(s) to HDF5.
    """

    with h5py.File(hdf5, "a") as h5:
        if "gui" not in h5.keys():
            h5.create_group("gui")

        group = h5["gui"]
        def store_dataset(name, dtype):
            """
            For storing datasets
            """
            if name not in kwargs:
                return

            if name not in group.keys():
                group.create_dataset(name, (), data=kwargs[name],
                                     dtype=dtype)
            else:
                group[name][()] = kwargs[name]


        store_dataset("input_rzplot_minr", "f8")
        store_dataset("input_rzplot_maxr", "f8")
        store_dataset("input_rzplot_numr", "i8")
        store_dataset("input_rzplot_minz", "f8")
        store_dataset("input_rzplot_maxz", "f8")
        store_dataset("input_rzplot_numz", "i8")


def retrieve(hdf5, **kwargs):
    """
    Retrieves parameter(s) from HDF5.

    kwargs : keys are names of the parameters to be retrieved and values
    are the default values to be used if the parameter is not found.
    """
    out = {}
    with h5py.File(hdf5, "r") as h5:
        if "gui" in h5.keys():
            for k in kwargs.keys():
                if k in h5["gui"]:
                    out[k] = h5["gui"][k][()]

    for k in kwargs.keys():
        if k not in out.keys():
            out[k] = kwargs[k]

    return out
