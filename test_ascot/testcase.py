"""
Test case module.
"""

import a5py.ascot5io.B_TC as B_TC
import a5py.ascot5io.E_TC as E_TC
import a5py.ascot5io.plasma_1D as plasma_1D
import a5py.ascot5io.wall_2D as wall_2D

def createbase(fn, Bxyz, Exyz, n, T):
    """
    Create test case base with uniform EM fields and plasma.

    fn : str
        Full path to HDF5 file.
    Bxyz : real 3 x 1 numpy array
        Magnetic field Bxyz.
    Exyz : real 3 x 1 numpy array
        Electric field Exyz.
    n : real
        Plasma density.
    T : real
        Plasma temperature (eV).
    """

    J = np.array([0,0,0,0,0,0,0,0,0])
    rhoval = 0.5
    B_TC.write_hdf5(fn, Bxyz, J, rhoval)
    
    E_TC.write_hdf5(fn, Exyz)

    Nrho = 2
    Nion = 1
    znum = 1
    anum = 1
    rho = np.array([0, 1])
    ndens = 0*ones(rho.shape)
    ntemp = 0*ones(rho.shape)
    edens = n*ones(rho.shape)
    etemp = T*ones(rho.shape)
    idens = n*ones(rho.shape)
    itemp = T*ones(rho.shape)
    plasma_1D.write_hdf5(fn, Nrho, Nion, znum, anum, rho, ndens, ntemp, edens, etemp, idens, itemp)

    n = 4
    R = [0.1, 100, 100, 0.1]
    z = [-100, -100, 100, 100] 
    wall2D.write_hdf5(fn, n, R, z)


def flagsToZero(options,flags):
    for i in options:
        if i.startswith(flags):
            setattr(options,i,0)
    return options


def compile(flags={}):
    """
    Compile ASCOT5.

    Parameters
    ----------

    flags : dictionary, optional
        Optional ascot5.h compiler options in dictionary 
        format [name, value].
    """
    
