import numpy as np
import h5py

def write_hdf5(fn, rlim, zlim, psirz, Br, Bphi, Bz, axisRz, psivals):
    """Write 2D magnetic field input in hdf5 file.

    Keyword arguments:
    fn      -- path to hdf5 file
    rlim    -- [r_min, r_max] values of the R grid
    zlim    -- [z_min, z_max] values of the z grid
    psirz   -- psi values in Rz grid
    Br      -- Br values in Rz grid
    Bphi    -- Bphi values in Rz grid
    Bz      -- Bz values in Rz grid
    axisRz  -- magnetic axis (R,z) coordinates
    psivals -- psi values at magnetic axis and X-point (separatrix)
    """
    
    # Create bfield group if one does not exists and set 
    # type to 2D. If one exists, only set type to 2D.
    f = h5py.File(fn, "a")
    if not "/bfield" in f:
        o = f.create_group('bfield')
        o.attrs["type"] = np.string_("2D")
    
    else:
        o = f["bfield"]
        del o.attrs["type"]
        o.attrs["type"] = np.string_("2D")
        
    # Remove 2D field if one is already present
    if  "/bfield/2D" in f:
        del f["/bfield/2D"]

    # Deduce field dimensions from psirz matrix
    n_r = psirz.shape[1];
    n_z = psirz.shape[0];

    f.create_group('bfield/2D')
    f.create_dataset('bfield/2D/r_min', data=rlim[0])
    f.create_dataset('bfield/2D/r_max', data=rlim[1])
    f.create_dataset('bfield/2D/n_r', data=n_r)

    f.create_dataset('bfield/2D/z_min', data=zlim[0])
    f.create_dataset('bfield/2D/z_max', data=zlim[1])
    f.create_dataset('bfield/2D/n_z', data=n_z)

    f.create_dataset('bfield/2D/psi', data=psirz.flatten(order='C'))
    f.create_dataset('bfield/2D/B_r', data=Br.flatten(order='C'))
    f.create_dataset('bfield/2D/B_phi', data=Bphi.flatten(order='C'))
    f.create_dataset('bfield/2D/B_z', data=Bz.flatten(order='C'))

    # Magnetic axis and psi values
    f.create_dataset('bfield/2D/axis_r', data=axisRz[0])
    f.create_dataset('bfield/2D/axis_z', data=axisRz[1])

    f.create_dataset('bfield/2D/psi0', data=psivals[0])
    f.create_dataset('bfield/2D/psi1', data=psivals[1])
    f.close()
