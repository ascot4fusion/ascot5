"""
Test case module.
"""
import numpy as np
import a5py.ascot5io.B_TC as B_TC
import a5py.ascot5io.E_TC as E_TC
import a5py.ascot5io.plasma_1D as plasma_1D
import a5py.ascot5io.wall_2D as wall_2D
import a5py.ascot5io.options as options

import sys
sys.path.append('..')
import opt

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
    znum = np.array([1])
    anum = np.array([1])
    rho = np.array([0, 1])
    ndens = 0*np.ones(rho.shape)
    ntemp = 0*np.ones(rho.shape)
    edens = n*np.ones(rho.shape)
    etemp = T*np.ones(rho.shape)
    idens = n*np.ones(rho.shape)
    itemp = T*np.ones(rho.shape)
    plasma_1D.write_hdf5(fn, Nrho, Nion, znum, anum, rho, ndens, ntemp, edens, etemp, idens, itemp)

    n = 4
    R = np.array([0.1, 100, 100, 0.1])
    z = np.array([-100, -100, 100, 100])
    wall_2D.write_hdf5(fn, n, R, z)

    o = opt.generateopt()
    o = flagsToZero(o,"ENABLE")
    o = flagsToZero(o,"ENDCOND")

    o["RECORD_GO_AS_GC"] = np.array([0],dtype='i4')

    options.write_hdf5(fn, o)

def flagsToZero(options,flags):
    for i in options:
        if i.startswith(flags):
            options[i] = np.array([0],dtype='i4')
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
    
