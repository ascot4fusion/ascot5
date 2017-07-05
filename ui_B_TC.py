import h5py
import numpy as np

def write_hdf5(fn, Bxyz, gradB):
	f = h5py.File(fn, "a")
	o = f.create_group("bfield")
        o.attrs["type"] = np.string_("B_TC");
	o = o.create_group("B_TC")
	o.create_dataset("Bxyz", (3,1), dtype='f8', data = Bxyz)
	o.create_dataset("gradB", (9,1), dtype='f8', data = gradB)
	f.close()

