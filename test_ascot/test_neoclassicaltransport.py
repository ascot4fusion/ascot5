"""
Test collisional transport in tokamak magnetic field.

This file initializes, runs, and plots test case for checking that ASCOT5
reproduces collsional transport in tokamak magnetic field, i.e., the
so-called neoclassical transport, correctly. This transport is diffusive
with the diffusion coefficient being determined by the ratio between
orbit bounce frequency and collision frequency. If collision frequency is a
lot higher than the bounce frequency, the transport is so-called
Pfirsch-Schlüter diffusion. If opposite is true, then the transport is called
banana transport. Between these two regimes is the plateau regime where
transport does not depend on collision frequency.

The neoclassical transport is tested on an analytical tokamak field that has
circular cross-section. Electrons are used as markers.

To init, run and check this test, call this script without any arguments. To
do only one of the above, call this script with an argument "init", "run", or
"check".

File: test_neoclassicaltransport.py
"""

import sys

import numpy                   as np
import scipy.constants         as constants
import matplotlib.pyplot       as plt

import a5py.ascot5io.ascot5    as ascot5
import a5py.ascot5io.options   as options
import a5py.ascot5io.B_GS      as B_GS
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.N0_3D     as N0_3D
import a5py.ascot5io.mrk_gc    as mrk

sys.path.insert(0, '../')
sys.path.insert(0, '.')
import opt
import test_ascot

from a5py.preprocessing.analyticequilibrium import psi0 as psifun

e       = constants.elementary_charge
m_e_AMU = constants.physical_constants["electron mass in u"][0]
m_e     = constants.physical_constants["electron mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

Nmrk  = 10
nscan = 3

Temp = 1e3
dens = np.power( 10, np.linspace(21.5, 23.5, nscan) )

R0 = 6.2
z0 = 0
Bphi0 = 5.3
psi_mult = 200
# ITER-like but circular equilibrium
psi_coeff = np.array([2.21808016e-02,  -1.28841781e-01,  -4.17718173e-02,
                      -6.22680280e-02,   6.20083978e-03,  -1.20524711e-03,
                      -3.70147050e-05,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      -0.155])
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
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-4
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_CLMBCOL_ENERGY"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10

    opt.settypes(odict)
    for i in range(1, nscan+1):
        options.write_hdf5(test_ascot.testfn, odict,
                           desc="NEOCLASS_GO" + str(i))

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
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-4
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_CLMBCOL_ENERGY"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8

    opt.settypes(odict)
    for i in range(1, nscan+1):
        options.write_hdf5(test_ascot.testfn, odict,
                           desc="NEOCLASS_GCF" + str(i))

    #**************************************************************************#
    #*                     Generate options for CLASS-GCA                      #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 5e-9
    odict["ADAPTIVE_TOL_COL"]          = 1e-1
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-4
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_CLMBCOL_ENERGY"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8

    opt.settypes(odict)
    for i in range(1, nscan+1):
        options.write_hdf5(test_ascot.testfn, odict,
                           desc="NEOCLASS_GCA" + str(i))

    #**************************************************************************#
    #*                    Marker input consisting of electrons                 #
    #*                                                                         #
    #**************************************************************************#
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = np.ones(ids.shape)
    mass   = m_e_AMU * np.ones(ids.shape)
    charge = 1       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    R      = 8       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    z      = 0       * np.ones(ids.shape)
    energy = 1e3     * np.ones(ids.shape)
    theta  = 2 * np.pi * np.random.rand(1,Nmrk)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    for i in range(1, nscan+1):
        mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time,
                       desc="NEOCLASS_GO" + str(i))
        mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time,
                       desc="NEOCLASS_GCF" + str(i))
        mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time,
                       desc="NEOCLASS_GCA" + str(i))

    #**************************************************************************#
    #*     Uniform magnetic field with values scanned from Bmin to Bmax        #
    #*                                                                         #
    #**************************************************************************#
    for i in range(1, nscan+1):
        B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                        desc="NEOCLASS_GO" + str(i))
        B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                        desc="NEOCLASS_GCF" + str(i))
        B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                        desc="NEOCLASS_GCA" + str(i))

    #**************************************************************************#
    #* Plasma consisting of protons only (to avoid e-e collisions)             #
    #*                                                                         #
    #**************************************************************************#
    Nrho  = 3
    Nion  = 1
    znum  = np.array([1])
    anum  = np.array([1])
    rho   = np.array([0, 0.5, 100])
    edens = 1    * np.ones(rho.shape)
    etemp = Temp * np.ones(rho.shape)
    itemp = 1e3  * np.ones(rho.shape)
    for i in range(1, nscan+1):
        idens = dens[i-1] * np.ones((rho.size, Nion))
        P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp,
                        desc="NEOCLASS_GO" + str(i))
        P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp,
                        desc="NEOCLASS_GCF" + str(i))
        P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp,
                        desc="NEOCLASS_GCA" + str(i))

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    for i in range(1, nscan+1):
        Exyz   = np.array([0, 0, 0])
        E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="NEOCLASS_GO" + str(i))
        E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="NEOCLASS_GCF" + str(i))
        E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="NEOCLASS_GCA" + str(i))

        nwall = 4
        Rwall = np.array([0.1, 100, 100, 0.1])
        zwall = np.array([-100, -100, 100, 100])
        W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                        desc="NEOCLASS_GO" + str(i))
        W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                        desc="NEOCLASS_GCF" + str(i))
        W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall,
                        desc="NEOCLASS_GCA" + str(i))

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
                         desc="NEOCLASS_GO" + str(i))
        N0_3D.write_hdf5(test_ascot.testfn,
                         N0Rmin, N0Rmax, N0nR,
                         N0zmin, N0zmax, N0nz,
                         N0pmin, N0pmax, N0np, N0dens,
                         desc="NEOCLASS_GCF" + str(i))
        N0_3D.write_hdf5(test_ascot.testfn,
                         N0Rmin, N0Rmax, N0nR,
                         N0zmin, N0zmax, N0nz,
                         N0pmin, N0pmax, N0np, N0dens,
                         desc="NEOCLASS_GCA" + str(i))

def run():
    """
    Run tests.
    """
    for test in ["NEOCLASS_GO", "NEOCLASS_GCF", "NEOCLASS_GCA"]:
        for i in range(1, nscan+1):
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

    R_omp    = np.linspace(6.2,8.2,1000)
    z_omp    = 0*np.ones(R_omp.shape)
    psi_omp  = psifun(R_omp/R0, z_omp/R0, psi_coeff[0], psi_coeff[1],
                      psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                      psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                      psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult
    psi_axis = psifun(R0/R0, 0, psi_coeff[0], psi_coeff[1],
                      psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                      psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                      psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult
    rho_omp  = np.sqrt(np.absolute( (psi_omp - psi_axis) / psi_axis ))

    DATA = {}
    for i in range(1, nscan+1):
        inistate = a5["NEOCLASS_GO" + str(i)].inistate.read()
        endstate = a5["NEOCLASS_GO" + str(i)].endstate.read()
        ri = np.interp(inistate["rho"], rho_omp, R_omp)
        rf = np.interp(endstate["rho"], rho_omp, R_omp)
        t  = endstate["time"]

        DATA["GO" + str(i)] = np.mean( np.power(ri - rf, 2) / t )

        inistate = a5["NEOCLASS_GCF" + str(i)].inistate.read()
        endstate = a5["NEOCLASS_GCF" + str(i)].endstate.read()
        ri = np.interp(inistate["rho"], rho_omp, R_omp)
        rf = np.interp(endstate["rho"], rho_omp, R_omp)
        t  = endstate["time"]

        DATA["GCF" + str(i)] = np.mean( np.power(ri - rf, 2) / t )

        inistate = a5["NEOCLASS_GCA" + str(i)].inistate.read()
        endstate = a5["NEOCLASS_GCA" + str(i)].endstate.read()
        ri = np.interp(inistate["rho"], rho_omp, R_omp)
        rf = np.interp(endstate["rho"], rho_omp, R_omp)
        t  = endstate["time"]

        DATA["GCA" + str(i)] = np.mean( np.power(ri - rf, 2) / t )

    DGO  = np.zeros( (nscan, 1) )
    DGCF = np.zeros( (nscan, 1) )
    DGCA = np.zeros( (nscan, 1) )
    Dps  = np.zeros( (nscan, 1) )
    Dp   = np.zeros( (nscan, 1) )
    Db   = np.zeros( (nscan, 1) )
    for i in range(1, nscan+1):
        DGO[i-1]  = DATA["GO" + str(i)]
        DGCF[i-1] = DATA["GCF" + str(i)]
        DGCA[i-1] = DATA["GCA" + str(i)]

        clog     = 14.7
        collfreq = 4.8e-14 * dens[i-1] * np.power(Temp,-3.0/2) * clog
        gyrolen  = np.sqrt(2 * m_e * 1e3 * e) / ( Bphi0 * e )
        eps      = np.power(1.8/6.2, 3.0/2.0)
        q        = 2.2
        v        = np.sqrt(2 * 1e3 * e / m_e)
        Dclass   = q * np.power(gyrolen, 2) * v / 6.2
        Dps[i-1] = collfreq * np.power(gyrolen,2)
        Dp[i-1]  = q * q * Dps[i-1] * 2 * 2 * 2 * 2
        Db[i-1]  = np.power(6.2/1.9, 3.0/2.0) * Dp[i-1]

    plt.figure()
    plt.plot(np.log10(dens), np.log10(DGO))
    plt.plot(np.log10(dens), np.log10(DGCF))
    plt.plot(np.log10(dens), np.log10(DGCA))
    plt.plot(np.log10(dens), np.log10(Dps))
    plt.plot(np.log10(dens), np.log10(Dp))
    plt.plot(np.log10(dens), np.log10(Db))
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
