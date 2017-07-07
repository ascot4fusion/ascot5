import h5py
import numpy as np

def write_hdf5(fn, Bxyz, gradB, axisr, axisz, psival, rhoval):
	f = h5py.File(fn, "a")
	o = f.create_group("bfield")
        o.attrs["type"] = np.string_("B_TC");
	o = o.create_group("B_TC")
        o.create_dataset("axisr", (1,), dtype='f8', data = axisr)
        o.create_dataset("axisz", (1,), dtype='f8', data = axisz)
        o.create_dataset("psival", (1,), dtype='f8', data = psival)
        o.create_dataset("rhoval", (1,), dtype='f8', data = rhoval)
	o.create_dataset("Bxyz", (3,1), dtype='f8', data = Bxyz)
	o.create_dataset("gradB", (9,1), dtype='f8', data = gradB)
	f.close()

