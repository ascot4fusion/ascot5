"""
Test collisional transport in uniform magnetic field.

This file initializes, runs, and plots test case for checking that ASCOT5
reproduces collsional transport in uniform magnetic field, i.e., the
so-called classical transport, correctly. This transport is diffusive
with the diffusion coefficient being inversely proportional to magnetic
field strength squared.

Tests are done in a Cartesian (B_TC) magnetic field using protons as markers.

To init, run and check this test, call this script without any arguments. To
do only one of the above, call this script with an argument "init", "run", or
"check".

File: test_classicaltransport.py
"""

import sys

import numpy                   as np
import scipy.constants         as constants
import matplotlib.pyplot       as plt

import a5py.ascot5io.ascot5    as ascot5
import a5py.ascot5io.options   as options
import a5py.ascot5io.B_TC      as B_TC
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.N0_3D     as N0_3D
import a5py.ascot5io.mrk_gc    as mrk

sys.path.insert(0, '../')
sys.path.insert(0, '.')
import opt
import test_ascot

e       = constants.elementary_charge
m_p_AMU = constants.physical_constants["proton mass in u"][0]
m_p     = constants.physical_constants["proton mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

Te   = 1e3
ne   = 1e22
Bmin = 1
Bmax = 10
nB   = 5

# Number of markers. Adjust to change the time it takes to simulate tests.
Nmrk = 100

def init():
    """
    Initialize tests

    This function initializes parameter scan for three test cases:
    - CLASS-GO tests gyro-orbit scheme
    - CLASS-GCF tests guiding center fixed-scheme
    - CLASS-GCA tests guiding center adaptive scheme

    Each test case is initialized nB times; each with an decreasing value for
    the magnetic field strength. The test cases are named with running index on
    their prefix e.g. CLASS-GO1, CLASS-GO2, and so on.

    All cases are run without energy collisions.

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                     Generate options for CLASS-GO                       #
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
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_CLMBCOL_ENERGY"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10

    opt.settypes(odict)
    for i in range(1, nB+1):
        options.write_hdf5(test_ascot.testfn, odict, desc="CLASS_GO" + str(i))

    #**************************************************************************#
    #*                     Generate options for CLASS-GCF                      #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-9
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_CLMBCOL_ENERGY"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8

    opt.settypes(odict)
    for i in range(1, nB+1):
        options.write_hdf5(test_ascot.testfn, odict, desc="CLASS_GCF" + str(i))

    #**************************************************************************#
    #*                     Generate options for CLASS-GCA                      #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 1e-9
    odict["ADAPTIVE_TOL_COL"]          = 1e-1
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_CLMBCOL_ENERGY"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8

    opt.settypes(odict)
    for i in range(1, nB+1):
        options.write_hdf5(test_ascot.testfn, odict, desc="CLASS_GCA" + str(i))

    #**************************************************************************#
    #*                    Marker input consisting of protons                   #
    #*                                                                         #
    #**************************************************************************#
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = np.ones(ids.shape)
    mass   = m_p_AMU * np.ones(ids.shape)
    charge = 1       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    R      = 5       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    z      = 0       * np.ones(ids.shape)
    energy = 1e3     * np.ones(ids.shape)
    theta  = 2 * np.pi * np.random.rand(1,Nmrk)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    for i in range(1, nB+1):
        mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time, desc="CLASS_GO" + str(i))
        mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time, desc="CLASS_GCF" + str(i))
        mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time, desc="CLASS_GCA" + str(i))

    #**************************************************************************#
    #*     Uniform magnetic field with values scanned from Bmin to Bmax        #
    #*                                                                         #
    #**************************************************************************#
    gradB  = np.array([0,0,0,0,0,0,0,0,0])
    rhoval = 1.5
    Bvals = np.linspace( 1 / ( Bmax * Bmax ), 1 / ( Bmin * Bmin ), nB )
    Bvals = 1 / np.sqrt(Bvals)
    for i in range(1, nB+1):
        Bxyz   = np.array([Bvals[i-1], 0, 0])
        B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval,
                        desc="CLASS_GO" + str(i))
        B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval,
                        desc="CLASS_GCF" + str(i))
        B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval,
                        desc="CLASS_GCA" + str(i))

    #**************************************************************************#
    #* Plasma consisting of electrons only (to avoid proton-proton collisions) #
    #*                                                                         #
    #**************************************************************************#
    Nrho  = 3
    Nion  = 1
    znum  = np.array([1])
    anum  = np.array([1])
    rho   = np.array([0, 0.5, 100])
    edens = ne  * np.ones(rho.shape)
    etemp = Te  * np.ones(rho.shape)
    idens = 1   * np.ones((rho.size, Nion))
    itemp = 1e3 * np.ones(rho.shape)
    for i in range(1, nB+1):
        P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp, desc="CLASS_GO" + str(i))
        P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp, desc="CLASS_GCF" + str(i))
        P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp, desc="CLASS_GCA" + str(i))

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    for i in range(1, nB+1):
        Exyz   = np.array([0, 0, 0])
        E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="CLASS_GO" + str(i))
        E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="CLASS_GCF" + str(i))
        E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="CLASS_GCA" + str(i))

        nwall = 4
        Rwall = np.array([0.1, 100, 100, 0.1])
        zwall = np.array([-100, -100, 100, 100])
        W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                        desc="CLASS_GO" + str(i))
        W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                        desc="CLASS_GCF" + str(i))
        W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                        desc="CLASS_GCA" + str(i))

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
                         N0pmin, N0pmax, N0np, N0dens,
                         desc="CLASS_GO" + str(i))
        N0_3D.write_hdf5(test_ascot.testfn,
                         N0Rmin, N0Rmax, N0nR,
                         N0zmin, N0zmax, N0nz,
                         N0pmin, N0pmax, N0np, N0dens,
                         desc="CLASS_GCF" + str(i))
        N0_3D.write_hdf5(test_ascot.testfn,
                         N0Rmin, N0Rmax, N0nR,
                         N0zmin, N0zmax, N0nz,
                         N0pmin, N0pmax, N0np, N0dens,
                         desc="CLASS_GCA" + str(i))


def run():
    """
    Run tests.
    """
    for test in ["CLASS_GO", "CLASS_GCF", "CLASS_GCA"]:
        for i in range(1, nB+1):
            test_ascot.set_and_run(test + str(i))


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

    DATA = {}
    for i in range(1, nB+1):
        inistate = a5["CLASS_GO" + str(i)].inistate.read()
        yi = inistate["R"] * np.sin( inistate["phi"] * np.pi / 180 )
        zi = inistate["z"]
        ti = inistate["time"]

        endstate = a5["CLASS_GO" + str(i)].endstate.read()
        ye = endstate["R"] * np.sin( endstate["phi"] * np.pi / 180 )
        ze = endstate["z"]
        te = endstate["time"]

        dr2 = ( yi - ye ) * ( yi - ye ) + ( zi - ze ) * ( zi - ze )
        DATA["GO" + str(i)] = np.mean( dr2 / ( te - ti ) )

        inistate = a5["CLASS_GCF" + str(i)].inistate.read()
        yi = inistate["R"] * np.sin( inistate["phi"] * np.pi / 180 )
        zi = inistate["z"]
        ti = inistate["time"]

        endstate = a5["CLASS_GCF" + str(i)].endstate.read()
        ye = endstate["R"] * np.sin( endstate["phi"] * np.pi / 180 )
        ze = endstate["z"]
        te = endstate["time"]

        dr2 = ( yi - ye ) * ( yi - ye ) + ( zi - ze ) * ( zi - ze )
        DATA["GCF" + str(i)] = np.mean( dr2 / ( te - ti ) )

        inistate = a5["CLASS_GCA" + str(i)].inistate.read()
        yi = inistate["R"] * np.sin( inistate["phi"] * np.pi / 180 )
        zi = inistate["z"]
        ti = inistate["time"]

        endstate = a5["CLASS_GCA" + str(i)].endstate.read()
        ye = endstate["R"] * np.sin( endstate["phi"] * np.pi / 180 )
        ze = endstate["z"]
        te = endstate["time"]

        dr2 = ( yi - ye ) * ( yi - ye ) + ( zi - ze ) * ( zi - ze )
        DATA["GCA" + str(i)] = np.mean( dr2 / ( te - ti ) )

    DGO  = np.zeros( (nB, 1) )
    DGCF = np.zeros( (nB, 1) )
    DGCA = np.zeros( (nB, 1) )
    for i in range(1, nB+1):
        DGO[i-1]  = DATA["GO" + str(i)]
        DGCF[i-1] = DATA["GCF" + str(i)]
        DGCA[i-1] = DATA["GCA" + str(i)]

    xvals    = np.linspace( 1 / ( Bmax * Bmax ), 1 / ( Bmin * Bmin ), nB )
    Bvals    = 1 / np.sqrt(xvals)

    clog     = 14.7
    collfreq = 4.8e-14 * ne * np.power( Te, -3.0 / 2 ) * clog
    gyrolen  = np.sqrt( 2 * m_p * 1e3 * e ) / ( e * Bvals )
    Dclass   = collfreq * gyrolen * gyrolen

    print(Dclass / DGO)

    plt.figure()
    plt.plot(xvals, DGO)
    plt.plot(xvals, DGCF)
    plt.plot(xvals, DGCA)
    plt.plot(xvals, Dclass)
    plt.show()


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

    if(   sys.argv[1] == "init" ):
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
