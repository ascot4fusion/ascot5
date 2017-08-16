import h5py
from pylab import *
import numpy as np

def write_hdf5(fn, x1x2x3, y1y2y3, z1z2z3, flag):
    f = h5py.File(fn, "a")
    if not "/wall" in f:
        o = f.create_group('wall')
        o.attrs["type"] = np.string_("3D")
    
    else:
        o = f["wall"]
        del o.attrs["type"]
        o.attrs["type"] = np.string_("3D")
    
    if  "/wall/3D" in f:
        del f["/wall/3D"]

    Nelements = size(flag)
        
    f.create_group('wall/3D')
    f.create_dataset('wall/3D/x1x2x3', (Nelements,3), dtype='f8', data=x1x2x3)
    f.create_dataset('wall/3D/y1y2y3', (Nelements,3), dtype='f8', data=y1y2y3)
    f.create_dataset('wall/3D/z1z2z3', (Nelements,3), dtype='f8', data=z1z2z3)
    f.create_dataset('wall/3D/flag', (Nelements,1), dtype='i4', data=flag)
    f['wall/3D'].attrs['n_elements'] = size(flag)
    f['wall/3D'].attrs['min_x'] = np.amin(x1x2x3)
    f['wall/3D'].attrs['max_x'] = np.amax(x1x2x3)
    f['wall/3D'].attrs['min_y'] = np.amin(y1y2y3)
    f['wall/3D'].attrs['max_y'] = np.amax(y1y2y3)
    f['wall/3D'].attrs['min_z'] = np.amin(z1z2z3)
    f['wall/3D'].attrs['max_z'] = np.amax(z1z2z3)

    f.close()
