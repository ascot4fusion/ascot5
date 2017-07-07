import h5py
import numpy as np

def write_hdf5_particles(fn, ids, anum, znum, r, phi, z, vr, vphi, vz, weight):
    f = h5py.File(fn, "a")
    if not "/markers" in f:
        f.create_group('markers')
    
    if  "/markers/particle" in f:
        del f["/markers/particle"]
        

    f.create_group('markers/particle')

    f['markers'].attrs['n_particle'] = np.size(ids)

    if not 'n_guiding_center' in f['markers'].attrs:
        f['markers'].attrs['n_guiding_center'] = 0

    if not 'n_field_line' in f['markers'].attrs:
        f['markers'].attrs['n_field_line'] = 0

    f.create_dataset('markers/particle/r', data=r)
    f.create_dataset('markers/particle/phi', data=phi)
    f.create_dataset('markers/particle/z', data=z)
    f.create_dataset('markers/particle/v_r', data=vr)
    f.create_dataset('markers/particle/v_phi', data=vphi)
    f.create_dataset('markers/particle/v_z', data=vz)
    f.create_dataset('markers/particle/anum', data=anum)
    f.create_dataset('markers/particle/znum', data=znum)
    f.create_dataset('markers/particle/weight', data=weight)
    f.create_dataset('markers/particle/id', data=ids)

    f.close()
