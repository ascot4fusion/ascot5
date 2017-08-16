import numpy as np
import h5py

def write_hdf5(fn, r, phi, z, br, bphi, bz, s, axis_r, axis_phi, axis_z, n_periods):
    
    # Create bfield group if one does not exists and set 
    # type to ST. If one exists, only set type to ST.
    f = h5py.File(fn, "a")
    if not "/bfield" in f:
        o = f.create_group('bfield')
        o.attrs["type"] = np.string_("B_ST")
    
    else:
        o = f["bfield"]
        del o.attrs["type"]
        o.attrs["type"] = np.string_("B_ST")
        
    # Remove stellarator field if one is already present
    if  "/bfield/B_ST" in f:
        del f["/bfield/B_ST"]
    
    f.create_group('bfield/B_ST')
    # Coordinates (TODO change these to just min_r, max_r, n_r etc.)
    f.create_dataset('bfield/B_ST/r', dtype='f8', data=r)
    f.create_dataset('bfield/B_ST/phi', dtype='f8', data=phi)
    f.create_dataset('bfield/B_ST/z', dtype='f8', data=z)

    # Magnetic field data
    f.create_dataset('bfield/B_ST/B_r', dtype='f8', data=br)
    f.create_dataset('bfield/B_ST/B_phi', dtype='f8', data=bphi)
    f.create_dataset('bfield/B_ST/B_z', dtype='f8', data=bz)
    f.create_dataset('bfield/B_ST/s', dtype='f8', data=s)

    # Magnetic axis (not implemented)
    f.create_dataset('bfield/B_ST/axis_r', dtype='f8', data=axis_r)
    f.create_dataset('bfield/B_ST/axis_phi', dtype='f8', data=axis_phi)
    f.create_dataset('bfield/B_ST/axis_z', dtype='f8', data=axis_z)

    # Toroidal periods
    f.create_dataset('bfield/B_ST/toroidalPeriods', dtype='i4', data=n_periods)
    
    f.close()
