import h5py
import numpy as np

def write_hdf5(fn, r, z):
    f = h5py.File(fn, "a")
    if not "/wall" in f:
        o = f.create_group('wall')
        o.attrs["type"] = np.string_("2D")
    
    else:
        o = f["wall"]
        del o.attrs["type"]
        o.attrs["type"] = np.string_("2D")
    
    if  "/wall/2D" in f:
        del f["/wall/2D"]

    f.create_group('wall/2D')
    f['wall/2D'].attrs['n'] = r.size
    f.create_dataset('wall/2D/r', data=r)
    f.create_dataset('wall/2D/z', data=z)

    f.close()
