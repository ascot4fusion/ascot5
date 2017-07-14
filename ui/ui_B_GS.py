import numpy as np
import h5py

def write_hdf5(fn, R0, z0, B_phi0, psi0, psi1, psi_mult, psi_coeff):
    """Write analytical tokamak magnetic field input in hdf5 file.

    Keyword arguments:
    fn        -- path to hdf5 file
    R0        -- axis R coordinate
    z0        -- axis z coordinate
    B_phi0    -- on-axis toroidal magnetic field
    psi0      -- psi at axis
    psi1      -- psi at separatrix
    psi_mult  -- psi multiplier
    psi_coeff -- numpy array of psi coefficients
    """

    # Create bfield group if one does not exists and set 
    # type to B_GS. If one exists, only set type to B_GS.
    f = h5py.File(fn, "a")
    if not "/bfield" in f:
        o = f.create_group('bfield')
        o.attrs["type"] = np.string_("B_GS")
    
    else:
        o = f["bfield"]
        del o.attrs["type"]
        o.attrs["type"] = np.string_("B_GS")
        
    # Remove B_GS field if one is already present
    if  "/bfield/B_GS" in f:
        del f["/bfield/B_GS"]

    f.create_group('bfield/B_GS')
    f.create_dataset('bfield/B_GS/R0', data=R0, dtype='f8')
    f.create_dataset('bfield/B_GS/z0', data=z0, dtype='f8')
    f.create_dataset('bfield/B_GS/B_phi0', data=B_phi0, dtype='f8')
    f.create_dataset('bfield/B_GS/psi0', data=psi0, dtype='f8')
    f.create_dataset('bfield/B_GS/psi1', data=psi1, dtype='f8')
    f.create_dataset('bfield/B_GS/psi_mult', data=psi_mult, dtype='f8')
    f.create_dataset('bfield/B_GS/psi_coeff', data=psi_coeff, dtype='f8')
