import h5py
import numpy as np
	
def write_hdf5(fn, Exyz):
    f = h5py.File(fn, "a")
    o = f.create_group("E_TC")
    o.attrs["type"] = np.string_("E_TC")
    o.create_dataset("Exyz", (3,1), dtype='f8', data = Exyz)
    f.close()
