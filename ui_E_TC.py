import h5py
import numpy as np
	
def write_hdf5(fn, Exyz):
	f = h5py.File(fn, "a")
	o = f.create_group("efield")
        o.attrs["type"] = np.string_("E_TC");
	o = o.create_group("TC")
	o.create_dataset("Exyz", (3,1), dtype='f16', data = Exyz)
	f.close()

write_hdf5("ascot.h5", np.array([1, 0, 0])) 
