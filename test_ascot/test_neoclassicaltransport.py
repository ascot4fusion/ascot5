"""
Test collisional transport in tokamak magnetic field.

This file initializes, runs, and plots test case for checking that ASCOT5
reproduces collisional transport in tokamak magnetic field, i.e., the
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
m_p     = constants.physical_constants["proton mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]
eps0    = constants.physical_constants["electric constant"][0]

Nmrk  = 200
nscan = 20

Ti   = 1e3
ni   = np.power( 10, np.linspace(18.5, 22.5, nscan) )
Ekin = 1e3

R0 = 6.2
r0 = 1.8
z0 = 0
B0 = 5.3
psi_mult = 200

# ITER-like but circular equilibrium
psi_coeff = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                       6.200e-03, -1.205e-03, -3.701e-05,  0,
                       0,          0,          0,          0,         -0.155])

def init():
    """
    Initialize tests

    This function initializes parameter scan for three test cases:
    - NEOCLASS-GO tests gyro-orbit scheme
    - NEOCLASS-GCF tests guiding center fixed-scheme
    - NEOCLASS-GCA tests guiding center adaptive scheme

    Each test case is initialized nscan times; each with a increasing value for
    the plasma density. The test cases are named with running index on their
    prefix e.g. NEOCLASS-GO1, NEOCLASS-GO2, and so on.

    All cases are run without energy collisions.

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                  Generate options for NEOCLASS-GO                       #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 3e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-4
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_ENERGY_CCOLL"]      = 1

    opt.settypes(odict)
    for i in range(1, nscan+1):
        options.write_hdf5(test_ascot.testfn, odict,
                           desc="NEOCLASS_GO" + str(i))

    #**************************************************************************#
    #*                  Generate options for NEOCLASS-GCF                      #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 5e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-4
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_ENERGY_CCOLL"]      = 1

    opt.settypes(odict)
    for i in range(1, nscan+1):
        options.write_hdf5(test_ascot.testfn, odict,
                           desc="NEOCLASS_GCF" + str(i))

    #**************************************************************************#
    #*                  Generate options for NEOCLASS-GCA                      #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 5e-9
    odict["ADAPTIVE_TOL_COL"]          = 9e-1
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-4
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["DISABLE_ENERGY_CCOLL"]      = 1

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
    R      = (R0+r0) * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    z      = 0       * np.ones(ids.shape)
    energy = Ekin    * np.ones(ids.shape)
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
    #*             Analytical ITER-like but circular magnetic field            #
    #*                                                                         #
    #**************************************************************************#
    for i in range(1, nscan+1):
        B_GS.write_hdf5(test_ascot.testfn, R0, z0, B0, psi_mult, psi_coeff,
                        desc="NEOCLASS_GO" + str(i))
        B_GS.write_hdf5(test_ascot.testfn, R0, z0, B0, psi_mult, psi_coeff,
                        desc="NEOCLASS_GCF" + str(i))
        B_GS.write_hdf5(test_ascot.testfn, R0, z0, B0, psi_mult, psi_coeff,
                        desc="NEOCLASS_GCA" + str(i))

    #**************************************************************************#
    #*     Plasma consisting of protons only (to avoid e-e collisions)         #
    #*                                                                         #
    #**************************************************************************#
    Nrho  = 3
    Nion  = 1
    znum  = np.array([1])
    anum  = np.array([1])
    rho   = np.array([0, 0.5, 100])
    edens = 1    * np.ones(rho.shape)
    etemp = Ti * np.ones(rho.shape)
    itemp = Ti  * np.ones(rho.shape)
    for i in range(1, nscan+1):
        idens = ni[i-1] * np.ones((rho.size, Nion))
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

    This function makes a one plot that shows numerically evaluated diffusion
    coefficients for all modes as a function of effective collision frequency
    (collision frequency divided by orbit transit frequency). For comparison,
    the analytical results for different regimes are also shown and also
    the boundaries of the regimes themselves.
    """

    a5 = ascot5.Ascot(test_ascot.testfn)

    # Map rho values to R outer mid-plane values
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

    #**************************************************************************#
    #*               Evaluate the numerical diffusion coefficients             #
    #*                                                                         #
    #**************************************************************************#
    DGO  = np.zeros(nscan)
    DGCF = np.zeros(nscan)
    DGCA = np.zeros(nscan)
    for i in range(1, nscan+1):
        inistate = a5["NEOCLASS_GO" + str(i)].inistate.read()
        endstate = a5["NEOCLASS_GO" + str(i)].endstate.read()
        ri = np.interp(inistate["rho"], rho_omp, R_omp)
        rf = np.interp(endstate["rho"], rho_omp, R_omp)
        t  = endstate["time"]

        DGO[i-1] = np.mean( np.power(ri - rf, 2) / (t) )

        inistate = a5["NEOCLASS_GCF" + str(i)].inistate.read()
        endstate = a5["NEOCLASS_GCF" + str(i)].endstate.read()
        ri = np.interp(inistate["rho"], rho_omp, R_omp)
        rf = np.interp(endstate["rho"], rho_omp, R_omp)
        t  = endstate["time"]

        DGCF[i-1] = np.mean( np.power(ri - rf, 2) / (t) )

        inistate = a5["NEOCLASS_GCA" + str(i)].inistate.read()
        endstate = a5["NEOCLASS_GCA" + str(i)].endstate.read()
        ri = np.interp(inistate["rho"], rho_omp, R_omp)
        rf = np.interp(endstate["rho"], rho_omp, R_omp)
        t  = endstate["time"]

        DGCA[i-1] = np.mean( np.power(ri - rf, 2) / (t) )

    #**************************************************************************#
    #*                  Evaluate the analytical estimate                       #
    #*                                                                         #
    #**************************************************************************#

    eps = r0 / R0
    q   = 2.2 # This is verified numerically
    B   = 5.3

    gamma  = 1 + Ekin * e / ( m_e * c * c )
    v      = np.sqrt(1.0 - 1.0 / ( gamma * gamma ) ) * c
    omegat = (v / (q * R0))
    rhog   = gamma * m_e * v / (B * e)

    clog     = 14.5
    density  = np.power( 10, np.linspace(np.log10(ni[0]) - 1,
                                         np.log10(ni[-1] + 1), 50) )
    collfreq = (np.sqrt(2/np.pi) / 3) \
               * np.power(e*e / ( 4*np.pi*eps0 ), 2) \
               * ( 4*np.pi / np.sqrt( m_e*np.power(Ti*e, 3) ) ) * density * clog
    veff = collfreq/omegat
    # Add values needed for plotting a continuous curve
    veff = np.append(veff, [1, np.power(eps, 3.0/2.0)])
    veff.sort()

    Dps = q * q * veff * omegat * np.power(rhog, 2) / 2
    Dp  = q * q * rhog * rhog * omegat \
          * np.ones(veff.shape) / 2
    Db  = np.power(eps, -3.0/2.0) * Dps

    # x coordinate for plotting the numerical coefficients
    collfreq = (np.sqrt(2/np.pi) / 3) \
               * np.power(e*e / ( 4*np.pi*eps0 ), 2) \
               * ( 4*np.pi / np.sqrt( m_e*np.power(Ti*e, 3) ) ) * ni * clog
    veff_x   = collfreq/omegat

    #**************************************************************************#
    #*                                  Plot                                   #
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

    h.plot(np.log10(np.power([eps, eps], 3.0/2.0)), np.array([-6, 0]), 'black')
    h.plot(np.log10(np.array([1, 1])), np.array([-6, 0]), 'grey')

    h.plot(np.log10(veff_x), np.log10(DGO), linestyle='none', marker='*',
           markersize=12, alpha=0.5)
    h.plot(np.log10(veff_x), np.log10(DGCF), linestyle='none', marker='.',
           markersize=10, alpha=0.5)
    h.plot(np.log10(veff_x), np.log10(DGCA), linestyle='none', marker='^',
           markersize=6, alpha=0.5)

    ind = veff >= 1
    h.plot(np.log10(veff[ind]), np.log10(Dps[ind]), 'black')
    ind = np.logical_and(veff >= np.power(eps, 3.0/2.0), veff <= 1)
    h.plot(np.log10(veff[ind]), np.log10(Dp[ind]), 'black')
    ind = veff <= np.power(eps, 3.0/2.0)
    h.plot(np.log10(veff[ind]), np.log10(Db[ind]), 'black')

    h.set_xlim(-2, 1)
    h.set_ylim(-5, -1)
    h.tick_params(axis='y', direction='out')
    h.tick_params(axis='x', direction='out')
    h.xaxis.set(ticks=np.linspace(-2, 1, 4))
    h.yaxis.set(ticks=np.linspace(-5, -1, 5),
                ticklabels=['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$',
                            '$10^{-1}$'])
    h.set(xlabel=r"$\log_{10}\nu^*$", ylabel=r"$D$ [m$^2$/s]")

    h.text(-0.98, -4.9, r"$\nu^*=\epsilon^{3/2}$", fontsize=10,
           bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
    h.text(-0.17, -4.9, r"$\nu^*=1$", fontsize=10,
           bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
    h.text(-1.5, -3.4, r"$D_{B}$", fontsize=10)
    h.text(-0.5, -3, r"$D_{P}$", fontsize=10)
    h.text(0.4, -2.1, r"$D_{PS}$", fontsize=10)

    plt.savefig("test_neoclassicaltransport.png", dpi=72)
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
