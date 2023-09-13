"""Store and retrieve GUI related parameters to/from HDF5
"""
import numpy as np
import h5py

class GUIparams():

    def __init__(self):
        self.params = {}

    def add(self, **kwargs):
        """Add parameters.
        """
        for k in kwargs.keys():
            self.params[k] = kwargs[k]

    def store(self, hdf5, params):
        """Stores parameter(s) to HDF5.

        Only stores the ones whose name is on the given list "params".
        """

        with h5py.File(hdf5, "a") as h5:
            if "gui" not in h5.keys():
                h5.create_group("gui")

            group = h5["gui"]
            def store_dataset(name, dtype):
                """For storing datasets
                """
                if name not in params:
                    return

                if name not in group.keys():
                    group.create_dataset(name, data=self.params[name].get(),
                                         dtype=dtype)
                else:
                    group[name][()] = self.params[name].get()

            S = h5py.string_dtype(encoding='utf-8', length=None)
            store_dataset("input_rzplot_minr", "f8")
            store_dataset("input_rzplot_maxr", "f8")
            store_dataset("input_rzplot_numr", "i8")
            store_dataset("input_rzplot_minz", "f8")
            store_dataset("input_rzplot_maxz", "f8")
            store_dataset("input_rzplot_numz", "i8")
            store_dataset("input_rzplot_qnt",  S)

    def retrieve(self, hdf5):
        """Retrieves parameter(s) from HDF5.
        """
        S = h5py.string_dtype(encoding='utf-8', length=None)
        try:
            with h5py.File(hdf5, "r") as h5:
                if "gui" in h5.keys():
                    for k in h5["gui"]:
                        if k in self.params.keys():
                            val = h5["gui"][k]
                            if val.dtype == S:
                                self.params[k].set(val.asstr()[()])
                            else:
                                self.params[k].set(val[()])

        except:
            pass
