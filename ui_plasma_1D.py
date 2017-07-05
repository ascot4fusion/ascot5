import h5py
import numpy as np

def write_hdf5(fn, znum, anum, rho, ndens, ntemp, edens, etemp, idens, itemp):
    f = h5py.File(fn, "a")

    f.create_group('plasma')

    f.create_dataset('plasma/Z_num', data=znum)
    f.create_dataset('plasma/A_mass', data=anum)
    f['plasma'].attrs['n_ions'] = idens.shape[1]
    f['plasma'].attrs['n_neutrals'] = ndens.shape[1]

    # 1D plasma properties
    f.create_dataset('plasma/1D/rho', data=rho)
    f.create_dataset('plasma/1D/temp_0', data=ntemp)
    f.create_dataset('plasma/1D/dens_0', data=ndens)
    f.create_dataset('plasma/1D/temp_e', data=etemp)
    f.create_dataset('plasma/1D/dens_e', data=edens)
    f.create_dataset('plasma/1D/temp_i', data=itemp)
    f.create_dataset('plasma/1D/dens_i', data=idens)
    f['plasma/1D'].attrs['n_rho'] = rho.size

    f.close();
