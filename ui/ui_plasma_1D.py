import h5py
import numpy as np

def write_hdf5(fn, Nrho, Nion, znum, anum, rho, ndens, ntemp, edens, etemp, idens, itemp):
    """Write 1D plasma input in hdf5 file.

    Keyword arguments:
    fn      -- path to hdf5 file
    Nrho    -- Number of rho grid points
    Nion    -- Number of ions
    znum    -- Ion charge number [Nion x 1]
    anum    -- Ion mass number [Nion x 1]
    rho     -- rho grid array
    ndens   -- neutral density (1/m^3) NOT IMPLEMNTED
    ntemp   -- neutral temperature (eV) NOT IMPLEMENTED
    edens   -- electron density (1/m^3) [Nrho x 1]
    etemp   -- electron temperature (eV) [Nrho x 1]
    idens   -- ion density (1/m^3) [Nrho x 1] 
    itemp   -- ion temperature (eV) [Nrho x Nion]
    """
    
    # neutrals are currently not implemented
    ndens = np.zeros((Nrho,1))
    ntemp = np.zeros((Nrho,1))
    Nneutral = 1

    # check that input is valid
    if anum.size != Nion or znum.size != Nion:
        raise Exception('Number of ions in input not consistent')

    if rho.size != Nrho or edens.size != Nrho or etemp.size != Nrho or itemp.size != Nrho:
        raise Exception('Number of rho grid points in input not consistent')

    if Nrho*Nion != idens.size:
        raise Exception('Ion density data is not consisten with Nrho and Nion')

    if etemp[0] < 1 or etemp[0] > 1e5 or itemp[0] < 1 or itemp[0] >1e5:
        print "Warning: Check that temperature is given in eV"

    # convert ion density matrix in 1D array (which is how it is stored in hdf5)
    if idens.size != (Nion,Nrho):
        idens.flatten("C")
    else:
        idens.flatten("F")

    f = h5py.File(fn, "a")

    # clear existing groups
    if not "/plasma" in f:
        o = f.create_group('plasma')
    else:
        o = f["plasma"]
        del o.attrs["type"]

    if  "/plasma/p1d" in f:
        del f["/plasma/P_1D"]

    o.attrs["type"] = np.string_("P_1D")
    o = o.create_group("P_1D")

    f.create_dataset('plasma/Z_num', (Nion,1), dtype='i4', data=znum)
    f.create_dataset('plasma/A_mass', (Nion,1), dtype='i4', data=anum)
    f['plasma'].attrs['n_ions'] = Nion
    f['plasma'].attrs['n_neutrals'] = Nneutral

    # 1D plasma properties
    f.create_dataset('plasma/P_1D/rho', (Nrho,1), dtype='f8', data=rho)
    f.create_dataset('plasma/P_1D/temp_0', (Nrho,1), dtype='f8', data=ntemp)
    f.create_dataset('plasma/P_1D/dens_0', (Nrho,1), dtype='f8', data=ndens)
    f.create_dataset('plasma/P_1D/temp_e', (Nrho,1), dtype='f8', data=etemp)
    f.create_dataset('plasma/P_1D/dens_e', (Nrho,1), dtype='f8', data=edens)
    f.create_dataset('plasma/P_1D/temp_i', (Nrho,1), dtype='f8', data=itemp)
    f.create_dataset('plasma/P_1D/dens_i', (Nrho*Nion,1), dtype='f8', data=idens)
    f['plasma/P_1D'].attrs['n_rho'] = Nrho

    f.close();
