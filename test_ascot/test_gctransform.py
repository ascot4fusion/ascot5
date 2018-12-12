"""
Test guiding center transformation.

File: test_gctransform.py
"""
import sys

import numpy                   as np
import scipy.constants         as constants
import matplotlib.pyplot       as plt

import a5py.ascot5io.ascot5    as ascot5
import a5py.ascot5io.orbits    as orbits
import a5py.ascot5io.options   as options
import a5py.ascot5io.B_GS      as B_GS
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.N0_3D     as N0_3D
import a5py.ascot5io.mrk_gc    as mrk

from a5py.preprocessing.analyticequilibrium import psi0 as psifun

sys.path.insert(0, '../')
sys.path.insert(0, '.')
import opt
import test_ascot

e       = constants.elementary_charge
m_e_AMU = constants.physical_constants["electron mass in u"][0]
m_e     = constants.physical_constants["electron mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

psi_mult  = 200
R0        = 6.2
z0        = 0
Bphi0     = 5.3
psi_coeff = np.array([ 8.629e-02,  3.279e-01,  5.268e-01, -2.366e-01,
                       3.825e-01, -3.573e-01, -1.484e-02,  1.506e-01,
                       7.428e-01, -4.447e-01, -1.084e-01,  1.281e-02, -0.155])

def init():
    """
    Initialize tests

    This function initializes three test cases:
    - ORBFOL-GO tests VPA algorithm used in integrating gyro-orbit motion
    - ORBFOL-GCF tests RK4 used in integrating guiding center motion with fixed
      time-step
    - ORBFOL-GCA tests Cash-Karp used in integrating guiding center motion with
      adaptive time-step

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                     Generate options for GC                             #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="GCTRANSFORM_GC")

    #**************************************************************************#
    #*                     Generate options for GO                             #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="GCTRANSFORM_GO")

    #**************************************************************************#
    #*                     Generate options for GO2GC                          #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["RECORD_GO_AS_GC"]           = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="GCTRANSFORM_GO2GC")

    #**************************************************************************#
    #*                 Marker input consisting of an electron                  #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 1
    ids    = 1
    weight = 1
    pitch  = 0.4
    mass   = m_e_AMU
    charge = -1
    time   = 0
    R      = 7.6
    phi    = 90
    z      = 0
    theta  = 2
    energy = 10e6
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GCTRANSFORM_GC")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GCTRANSFORM_GO")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GCTRANSFORM_GO2GC")

    #**************************************************************************#
    #*                     Construct ITER-like magnetic field                  #
    #*                                                                         #
    #**************************************************************************#
    Rmin = 4; Rmax = 8.5; nR = 120; zmin = -4; zmax = 4; nz = 200;
    B_GS.write_hdf5_B_2D(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_GC")
    B_GS.write_hdf5_B_2D(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_GO")
    B_GS.write_hdf5_B_2D(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_GO2GC")

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="GCTRANSFORM_GC")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="GCTRANSFORM_GO")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="GCTRANSFORM_GO2GC")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_GC")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_GO")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_GO2GC")

    N0Rmin = 0
    N0Rmax = 100
    N0nR   = 2
    N0zmin = -100
    N0zmax = 100
    N0nz   = 2
    N0pmin = 0
    N0pmax = 2*np.pi
    N0np   = 2
    N0dens = np.array([ [ [0,0] , [0,0] ], [ [0,0] , [0,0] ] ])
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="GCTRANSFORM_GC")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="GCTRANSFORM_GO")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="GCTRANSFORM_GO2GC")

    Nrho  = 3
    Nion  = 1
    znum  = np.array([1])
    anum  = np.array([1])
    rho   = np.array([0, 0.5, 100])
    edens = 1e20 * np.ones(rho.shape)
    etemp = 1e3  * np.ones(rho.shape)
    idens = 1e20 * np.ones((rho.size, Nion))
    itemp = 1e3  * np.ones(rho.shape)
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_GC")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_GO")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_GO2GC")

def run():
    """
    Run tests.
    """
    for test in ["GCTRANSFORM_GC", "GCTRANSFORM_GO", "GCTRANSFORM_GO2GC"]:
        test_ascot.set_and_run(test)

def check():
    """
    Plot the results of these tests.

    This function makes four plots.
    - One that shows conservation of energy for all cases
    - One that shows conservation of magnetic moment for all cases
    - One that shows conservation of toroidal canonical momentum for all cases
    - And one that shows trajectories on a Rz plane for all cases
    """
    a5 = ascot5.Ascot(test_ascot.testfn)

    GCTRAN ={}
    go    = a5["GCTRANSFORM_GO"]["fo"]
    go2gc = a5["GCTRANSFORM_GO2GC"]["gc"]
    gc    = a5["GCTRANSFORM_GC"]["gc"]
    govpar = go["vpar"]
    gomu   = go["mu"]

    if True:
        plt.figure()
        plt.plot(go["R"], go["z"])
        plt.plot(go2gc["R"], go2gc["z"])
        plt.plot(gc["R"], gc["z"])
    else:
        plt.figure()
        x,y = pol2cart(go["R"], go["phi"]*np.pi/180)
        plt.plot(x,go["z"])
        x,y = pol2cart(go2gc["R"], go2gc["phi"]*np.pi/180)
        plt.plot(x,go2gc["z"])
        x,y = pol2cart(gc["R"], gc["phi"]*np.pi/180)
        plt.plot(x,gc["z"])
        plt.axis('equal')

    plt.figure()
    plt.plot(go["time"], gomu)
    plt.plot(go2gc["time"], go2gc["mu"])
    plt.plot(gc["time"], gc["mu"])

    plt.figure()
    plt.plot(go["time"], govpar)
    plt.plot(go2gc["time"], go2gc["vpar"])
    plt.plot(gc["time"], gc["vpar"])

    plt.show()

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

if __name__ == '__main__':
    if( len(sys.argv) == 1 ):
        print("Initializing tests.")
        init()
        print("Initialization complete.")
        print("")
        print("Running tests.")
        run()
        print("Runs complete.")
        print("")
        print("Checking test results.")
        check()
        print("Testing complete.")
        sys.exit()

    if(len(sys.argv) > 2):
        print("Too many arguments.")
        print("Only \"init\", \"run\" or \"check\" is accepted.")
        print("Aborting.")
        sys.exit()

    if( sys.argv[1] == "init" ):
        print("Initializing tests.")
        init()
        print("Initialization complete.")
        sys.exit()

    elif( sys.argv[1] == "run" ):
        print("Running tests.")
        run()
        print("Runs complete.")
        sys.exit()

    elif( sys.argv[1] == "check" ):
        print("Checking test results.")
        check()
        print("Testing complete.")
        sys.exit()

    else:
        print("Too many arguments.")
        print("Only \"init\", \"run\" or \"check\" is accepted.")
        print("Aborting.")
        sys.exit()
