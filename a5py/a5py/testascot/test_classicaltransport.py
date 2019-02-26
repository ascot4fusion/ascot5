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
import a5py.ascot5io.mrk_gc    as mrk

import a5py.testascot.helpers as helpers

e       = constants.elementary_charge
m_p_AMU = constants.physical_constants["proton mass in u"][0]
m_p     = constants.physical_constants["proton mass"][0]
m_e     = constants.physical_constants["electron mass"][0]
eps0    = constants.physical_constants["electric constant"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

Te   = 1e3
ne   = 1e22
Bmin = 1
Bmax = 10
nB   = 6

# Number of markers. Adjust to change the time it takes to simulate tests.
Nmrk = 200

def init():
    """
    Initialize tests

    This function initializes parameter scan for three test cases:
    - CLASS_GO tests gyro-orbit scheme
    - CLASS_GCF tests guiding center fixed-scheme
    - CLASS_GCA tests guiding center adaptive scheme

    Each test case is initialized nB times; each with an decreasing value for
    the magnetic field strength. The test cases are named with running index on
    their prefix e.g. CLASS_GO1, CLASS_GO2, and so on.

    All cases are run without energy collisions.

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                     Generate options for CLASS-GO                       #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1

    for i in range(1, nB+1):
        options.write_hdf5(helpers.testfn, odict, desc="CLASS_GO" + str(i))

    #**************************************************************************#
    #*                     Generate options for CLASS-GCF                      #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-9
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1

    for i in range(1, nB+1):
        options.write_hdf5(helpers.testfn, odict, desc="CLASS_GCF" + str(i))

    #**************************************************************************#
    #*                     Generate options for CLASS-GCA                      #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 1e-8
    odict["ADAPTIVE_TOL_COL"]          = 1e-1
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1

    for i in range(1, nB+1):
        options.write_hdf5(helpers.testfn, odict, desc="CLASS_GCA" + str(i))

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
    phi    = 0       * np.ones(ids.shape)
    z      = 0       * np.ones(ids.shape)
    energy = 1e3     * np.ones(ids.shape)
    theta  = 2 * np.pi * np.random.rand(1,Nmrk)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    for i in range(1, nB+1):
        mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time, desc="CLASS_GO" + str(i))
        mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, energy, pitch, theta,
                       weight, time, desc="CLASS_GCF" + str(i))
        mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
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
        B_TC.write_hdf5(helpers.testfn, Bxyz, gradB, rhoval,
                        desc="CLASS_GO" + str(i))
        B_TC.write_hdf5(helpers.testfn, Bxyz, gradB, rhoval,
                        desc="CLASS_GCF" + str(i))
        B_TC.write_hdf5(helpers.testfn, Bxyz, gradB, rhoval,
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
        P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp, desc="CLASS_GO" + str(i))
        P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp, desc="CLASS_GCF" + str(i))
        P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                        edens, etemp, idens, itemp, desc="CLASS_GCA" + str(i))

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    for i in range(1, nB+1):
        Exyz   = np.array([0, 0, 0])
        E_TC.write_hdf5(helpers.testfn, Exyz, desc="CLASS_GO" + str(i))
        E_TC.write_hdf5(helpers.testfn, Exyz, desc="CLASS_GCF" + str(i))
        E_TC.write_hdf5(helpers.testfn, Exyz, desc="CLASS_GCA" + str(i))

        nwall = 4
        Rwall = np.array([0.1, 100, 100, 0.1])
        zwall = np.array([-100, -100, 100, 100])
        W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                        desc="CLASS_GO" + str(i))
        W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                        desc="CLASS_GCF" + str(i))
        W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                        desc="CLASS_GCA" + str(i))

        helpers.write_N0_3D_dummy(helpers.testfn, desc="CLASS_GO" + str(i))
        helpers.write_N0_3D_dummy(helpers.testfn, desc="CLASS_GCF" + str(i))
        helpers.write_N0_3D_dummy(helpers.testfn, desc="CLASS_GCA" + str(i))


def run():
    """
    Run tests.
    """
    for test in ["CLASS_GO", "CLASS_GCF", "CLASS_GCA"]:
        for i in range(1, nB+1):
            helpers.set_and_run(test + str(i))


def check():
    """
    Plot the results of these tests.

    This function makes one plot that shows how diffusion coefficient evaluated
    numerically with different modes scales with magnetic field strength, and
    compares that to the classical diffusion coefficient.
    """
    a5 = ascot5.Ascot(helpers.testfn)

    # Diffusion occurs on a 2D plane
    ndim = 2

    # Evaluate diffusion coefficient according to equation
    # D = var(Delta x) / (2 * ndim * Delta t)
    # Because E[Delta x] = 0, var[Delta x] = E[(Delta x)^2]
    DGO  = np.zeros(nB)
    DGCF = np.zeros(nB)
    DGCA = np.zeros(nB)
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
        DGO[i-1] = np.mean( dr2 / ( te - ti ) ) / (2*ndim)

        inistate = a5["CLASS_GCF" + str(i)].inistate.read()
        yi = inistate["R"] * np.sin( inistate["phi"] * np.pi / 180 )
        zi = inistate["z"]
        ti = inistate["time"]

        endstate = a5["CLASS_GCF" + str(i)].endstate.read()
        ye = endstate["R"] * np.sin( endstate["phi"] * np.pi / 180 )
        ze = endstate["z"]
        te = endstate["time"]

        dr2 = ( yi - ye ) * ( yi - ye ) + ( zi - ze ) * ( zi - ze )
        DGCF[i-1] = np.mean( dr2 / ( te - ti ) ) / (2*ndim)

        inistate = a5["CLASS_GCA" + str(i)].inistate.read()
        yi = inistate["R"] * np.sin( inistate["phi"] * np.pi / 180 )
        zi = inistate["z"]
        ti = inistate["time"]

        endstate = a5["CLASS_GCA" + str(i)].endstate.read()
        ye = endstate["R"] * np.sin( endstate["phi"] * np.pi / 180 )
        ze = endstate["z"]
        te = endstate["time"]

        dr2 = ( yi - ye ) * ( yi - ye ) + ( zi - ze ) * ( zi - ze )
        DGCA[i-1] = np.mean( dr2 / ( te - ti ) ) / (2*ndim)

    #**************************************************************************#
    #*                  Evaluate the analytical estimate                       #
    #*                                                                         #
    #**************************************************************************#

    clog     = 13.4
    collfreq = (m_e/m_p) * (np.sqrt(2/np.pi)/3) \
               * np.power(e*e/(4*np.pi*eps0), 2) \
               * (4*np.pi / np.sqrt(m_e*np.power(Te*e,3) ) ) * ne * clog

    xvals    = np.linspace( 1 / ( Bmax * Bmax ), 1 / ( Bmin * Bmin ), nB )
    Bvals    = 1 / np.sqrt(xvals)
    gamma  = 1 + 1e3*e / ( m_p * c * c )
    v      = np.sqrt(1.0 - 1.0 / ( gamma * gamma ) ) * c
    rhog   = gamma * m_p * v / (Bvals * e)
    Dclass = collfreq * rhog * rhog / 2

    #**************************************************************************#
    #*                            Plot                                         #
    #*                                                                         #
    #**************************************************************************#

    f = plt.figure(figsize=(11.9/2.54, 5/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    h = plt.gca()
    h.set_position([0.15, 0.25, 0.82, 0.7], which='both')

    scaling = 100
    h.plot(xvals, Dclass * scaling, 'black', linewidth=1, linestyle='-',
           label='Analytical')
    h.plot(xvals, DGO * scaling, linestyle='none', marker='*', markersize=11,
           alpha=0.5, label='GO')
    h.plot(xvals, DGCF * scaling, linestyle='none', marker='.', markersize=10,
           alpha=0.5, label='GCF')
    h.plot(xvals, DGCA * scaling, linestyle='none', marker='^', markersize=6,
           alpha=0.5, label='GCA')

    h.tick_params(axis='y', direction='out')
    h.tick_params(axis='x', direction='out')
    h.xaxis.set(ticks=np.linspace(0,1,11),
                ticklabels=[0, '', '', '', '', 5, '', '', '', '', 10])
    h.yaxis.set(ticks=np.linspace(0,8,9),
                ticklabels=[0, '', 2, '',  4, '', 6, '', 8])
    h.set(xlabel="$1/B^2$ [T$^{-2}$]", ylabel="$D$ [m$^2$/s]")

    h.legend(loc='upper left', frameon=False, fontsize=10, numpoints=1)

    plt.savefig("test_classicaltransport.png", dpi=300)
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
