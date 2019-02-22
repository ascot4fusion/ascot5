"""
Test conservation properties of the integration of the Hamiltonian motion.

This file initializes, runs, and plots test case for checking that ASCOT5
orbit integrators have the following properties:

1. gyro-orbit integrator (VPA) conserves energy, and magnetic moment and
   canonical toroidal momentum are bounded
2. guiding center fixed step scheme (RK4) conserves magnetic moment
   (trivially) and drift in energy and canonical toroidal momentum decrese
   with time-step
3. guiding center adaptive step scheme (Cash-Karp) conserves magnetic moment
   (trivially) and drift in energy and canonical toroidal momentum decrese
   with error-tolerance

Tests are done in ITER-like but analytic (B_GS) magnetic field without
collisions. Test particle is a energetic electron (and positron) so tests also
verify that ASCOT5 is valid in relativistic regime.

To init, run and check this test, call this script without any arguments. To
do only one of the above, call this script with an argument "init", "run", or
"check".

File: test_orbitfollowing.py
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

import a5py.testascot.helpers as helpers

from a5py.preprocessing.analyticequilibrium import psi0 as psifun

e       = constants.elementary_charge
m_e_AMU = constants.physical_constants["electron mass in u"][0]
m_e     = constants.physical_constants["electron mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

psi_mult  = 200
R0        = 6.2
z0        = 0
Bphi0     = 5.3

# ITER-like but circular equilibrium
psi_coeff = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                       6.200e-03, -1.205e-03, -3.701e-05,  0,
                       0,          0,          0,          0,         -0.155])

# Change this if you wish to use spline-interpolated magnetic field instead
use_spline = True

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
    #*                     Generate options for ORBFOL-GO                      #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-11
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10
    odict["ORBITWRITE_MAXPOINTS"]      = 50002

    options.write_hdf5(helpers.testfn, odict, desc="ORBFOL_GO")

    #**************************************************************************#
    #*                     Generate options for ORBFOL-GCF                     #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8
    odict["ORBITWRITE_MAXPOINTS"]      = 502

    options.write_hdf5(helpers.testfn, odict, desc="ORBFOL_GCF")

    #**************************************************************************#
    #*                     Generate options for ORBFOL-GCA                     #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 1e-10
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8
    odict["ORBITWRITE_MAXPOINTS"]      = 502

    options.write_hdf5(helpers.testfn, odict, desc="ORBFOL_GCA")

    #**************************************************************************#
    #*           Marker input consisting of an electron and positron           #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 2
    ids    = np.array([1, 2])
    weight = np.array([1, 1])
    pitch  = np.array([0.4, 0.9])
    mass   = m_e_AMU * np.array([1, 1])
    charge = 1       * np.array([1,-1])
    time   = 0       * np.array([1, 1])
    R      = 7.6     * np.array([1, 1])
    phi    = 90      * np.array([1, 1])
    z      = 0       * np.array([1, 1])
    theta  = 2       * np.array([1, 1])
    energy = 10e6    * np.array([1, 1])
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="ORBFOL_GO")
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="ORBFOL_GCF")
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="ORBFOL_GCA")

    #**************************************************************************#
    #*                     Construct ITER-like magnetic field                  #
    #*                                                                         #
    #**************************************************************************#
    if use_spline:
        B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                        desc="ORBFOL_GO")
        B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                        desc="ORBFOL_GCF")
        B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                        desc="ORBFOL_GCA")
    else:
        Rmin = 4; Rmax = 8.5; nR = 120; zmin = -4; zmax = 4; nz = 200;
        B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult,
                             psi_coeff, Rmin, Rmax, nR, zmin, zmax, nz,
                             desc="ORBFOL_GO")
        B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult,
                             psi_coeff, Rmin, Rmax, nR, zmin, zmax, nz,
                             desc="ORBFOL_GCF")
        B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult,
                             psi_coeff, Rmin, Rmax, nR, zmin, zmax, nz,
                             desc="ORBFOL_GCA")

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="ORBFOL_GO")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="ORBFOL_GCF")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="ORBFOL_GCA")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="ORBFOL_GO")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="ORBFOL_GCF")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="ORBFOL_GCA")

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
    N0_3D.write_hdf5(helpers.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="ORBFOL_GO")
    N0_3D.write_hdf5(helpers.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="ORBFOL_GCF")
    N0_3D.write_hdf5(helpers.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="ORBFOL_GCA")

    Nrho  = 3
    Nion  = 1
    znum  = np.array([1])
    anum  = np.array([1])
    rho   = np.array([0, 0.5, 100])
    edens = 1e20 * np.ones(rho.shape)
    etemp = 1e3  * np.ones(rho.shape)
    idens = 1e20 * np.ones((rho.size, Nion))
    itemp = 1e3  * np.ones(rho.shape)
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="ORBFOL_GO")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="ORBFOL_GCF")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="ORBFOL_GCA")


def run():
    """
    Run tests.
    """
    for test in ["ORBFOL_GO", "ORBFOL_GCF", "ORBFOL_GCA"]:
        helpers.set_and_run(test)

def check():
    """
    Plot the results of these tests.

    This function makes four plots.
    - One that shows conservation of energy for all cases
    - One that shows conservation of magnetic moment for all cases
    - One that shows conservation of toroidal canonical momentum for all cases
    - And one that shows trajectories on a Rz plane for all cases
    """
    a5 = ascot5.Ascot(helpers.testfn)

    f = plt.figure(figsize=(11.9/2.54, 8/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    h1 = f.add_subplot(1,4,1)
    h1.set_position([0.12, 0.72, 0.4, 0.25], which='both')

    h2 = f.add_subplot(1,4,2)
    h2.set_position([0.12, 0.44, 0.4, 0.25], which='both')

    h3 = f.add_subplot(1,4,3)
    h3.set_position([0.12, 0.155, 0.4, 0.25], which='both')

    h4 = f.add_subplot(1,4,4)
    h4.set_position([0.6, 0.3, 0.45, 0.45], which='both')

    #colors = ["#373e02", "#0a481e", "#03719c", "#0165fc", "#7e1e9c", "#cea2fd"]
    colors = ["b", "dodgerblue", "darkgreen", "forestgreen", "r", "tomato"]

    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for ORBFOL_GO             #
    #*                                                                         #
    #**************************************************************************#
    ORBFOL = {}
    ORBFOL["GO"] = {}
    orb = a5["ORBFOL_GO"]["orbits"].read()

    B = np.sqrt( orb["B_R"] * orb["B_R"] + orb["B_phi"] * orb["B_phi"] +
                 orb["B_z"] * orb["B_z"] )

    psi = psifun(orb["R"]/R0, orb["z"]/R0, psi_coeff[0], psi_coeff[1],
                 psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                 psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                 psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult

    vnorm = np.sqrt( orb["v_R"]   * orb["v_R"] +
                     orb["v_phi"] * orb["v_phi"] +
                     orb["v_z"]   * orb["v_z"] )

    vpar  = ( orb["v_R"] * orb["B_R"] + orb["v_phi"] * orb["B_phi"] +
              orb["v_z"] * orb["B_z"] ) / B

    gamma = np.sqrt(1 / ( 1 - vnorm * vnorm / (c * c) ) )

    ORBFOL["GO"]["time"] = orb["time"]
    ORBFOL["GO"]["id"]   = orb["id"]
    ORBFOL["GO"]["R"]    = orb["R"]
    ORBFOL["GO"]["z"]    = orb["z"]
    ORBFOL["GO"]["Ekin"] = (gamma - 1) * m_e * c * c
    ORBFOL["GO"]["mu"]   = ( ( m_e * gamma * gamma ) / ( 2 * B ) ) * \
                           ( vnorm * vnorm - vpar * vpar )
    ORBFOL["GO"]["ctor"] = gamma * m_e * orb["R"] * orb["v_phi"] + \
                           orb["charge"] * e * psi

    id1 = ORBFOL["GO"]["id"] == 1
    id2 = ORBFOL["GO"]["id"] == 2
    plot_relerr(h1, ORBFOL["GO"]["time"][id1], ORBFOL["GO"]["Ekin"][id1],
                colors[0])
    plot_relerr(h1, ORBFOL["GO"]["time"][id2], ORBFOL["GO"]["Ekin"][id2],
                colors[1])
    plot_relerr(h2, ORBFOL["GO"]["time"][id2], ORBFOL["GO"]["mu"][id2],
                colors[0])
    plot_relerr(h2, ORBFOL["GO"]["time"][id1], ORBFOL["GO"]["mu"][id1],
                colors[1])
    plot_relerr(h3, ORBFOL["GO"]["time"][id1], ORBFOL["GO"]["ctor"][id1],
                colors[0])
    plot_relerr(h3, ORBFOL["GO"]["time"][id2], ORBFOL["GO"]["ctor"][id2],
                colors[1])
    h4.plot(        ORBFOL["GO"]["R"][id1],    ORBFOL["GO"]["z"][id1],
                    colors[0])
    h4.plot(        ORBFOL["GO"]["R"][id2],    ORBFOL["GO"]["z"][id2],
                    colors[1])

    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for ORBFOL_GCF            #
    #*                                                                         #
    #**************************************************************************#
    ORBFOL["GCF"] = {}
    orb = a5["ORBFOL_GCF"]["orbits"].read()

    B = np.sqrt(np.power(orb["B_R"],2) + np.power(orb["B_phi"],2) +
                np.power(orb["B_z"],2))

    psi = psifun(orb["R"]/R0, orb["z"]/R0, psi_coeff[0], psi_coeff[1],
                 psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                 psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                 psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult

    # Note that mu is in eV / T
    gamma = np.sqrt( ( 1 + 2 * orb["mu"] * e * B / ( m_e * c * c ) ) /
                     ( 1 - orb["vpar"] * orb["vpar"] / ( c * c ) ) )

    ORBFOL["GCF"]["time"] = orb["time"]
    ORBFOL["GCF"]["id"]   = orb["id"]
    ORBFOL["GCF"]["R"]    = orb["R"]
    ORBFOL["GCF"]["z"]    = orb["z"]
    ORBFOL["GCF"]["Ekin"] = (gamma - 1) * m_e * c * c
    ORBFOL["GCF"]["mu"]   = orb["mu"] * e
    ORBFOL["GCF"]["ctor"] = gamma * m_e * orb["R"] * orb["vpar"] + \
                            orb["charge"] * e * psi

    id1 = ORBFOL["GCF"]["id"] == 1
    id2 = ORBFOL["GCF"]["id"] == 2
    plot_relerr(h1, ORBFOL["GCF"]["time"][id1], ORBFOL["GCF"]["Ekin"][id1],
                colors[2])
    plot_relerr(h1, ORBFOL["GCF"]["time"][id2], ORBFOL["GCF"]["Ekin"][id2],
                colors[3])
    plot_relerr(h2, ORBFOL["GCF"]["time"][id1], ORBFOL["GCF"]["mu"][id1],
                colors[2])
    plot_relerr(h2, ORBFOL["GCF"]["time"][id2], ORBFOL["GCF"]["mu"][id2],
                colors[3])
    plot_relerr(h3, ORBFOL["GCF"]["time"][id1], ORBFOL["GCF"]["ctor"][id1],
                colors[2])
    plot_relerr(h3, ORBFOL["GCF"]["time"][id2], ORBFOL["GCF"]["ctor"][id2],
                colors[3])
    h4.plot(        ORBFOL["GCF"]["R"][id1],    ORBFOL["GCF"]["z"][id1],
                colors[2])
    h4.plot(        ORBFOL["GCF"]["R"][id2],    ORBFOL["GCF"]["z"][id2],
                colors[3])

    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for ORBFOL_GCA            #
    #*                                                                         #
    #**************************************************************************#
    ORBFOL["GCA"] = {}
    orb = a5["ORBFOL_GCA"]["orbits"].read()

    B = np.sqrt(np.power(orb["B_R"],2) + np.power(orb["B_phi"],2) +
                np.power(orb["B_z"],2))

    psi = psifun(orb["R"]/R0, orb["z"]/R0, psi_coeff[0], psi_coeff[1],
                 psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                 psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                 psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult

    # Note that mu is in eV / T
    gamma = np.sqrt( ( 1 + 2 * orb["mu"] * e * B / ( m_e * c * c ) ) /
                     ( 1 - orb["vpar"] * orb["vpar"] / ( c * c ) ) )

    ORBFOL["GCA"]["time"] = orb["time"]
    ORBFOL["GCA"]["id"]   = orb["id"]
    ORBFOL["GCA"]["R"]    = orb["R"]
    ORBFOL["GCA"]["z"]    = orb["z"]
    ORBFOL["GCA"]["Ekin"] = (gamma - 1) * m_e * c * c
    ORBFOL["GCA"]["mu"]   = orb["mu"] * e
    ORBFOL["GCA"]["ctor"] = gamma * m_e * orb["R"] * orb["vpar"] + \
                            orb["charge"] * e * psi

    id1 = ORBFOL["GCA"]["id"] == 1
    id2 = ORBFOL["GCA"]["id"] == 2
    plot_relerr(h1, ORBFOL["GCA"]["time"][id1], ORBFOL["GCA"]["Ekin"][id1],
                colors[4])
    plot_relerr(h1, ORBFOL["GCA"]["time"][id2], ORBFOL["GCA"]["Ekin"][id2],
                colors[5])
    plot_relerr(h2, ORBFOL["GCA"]["time"][id1], ORBFOL["GCA"]["mu"][id1],
                colors[4])
    plot_relerr(h2, ORBFOL["GCA"]["time"][id2], ORBFOL["GCA"]["mu"][id2],
                colors[5])
    plot_relerr(h3, ORBFOL["GCA"]["time"][id1], ORBFOL["GCA"]["ctor"][id1],
                colors[4])
    plot_relerr(h3, ORBFOL["GCA"]["time"][id2], ORBFOL["GCA"]["ctor"][id2],
                colors[5])
    h4.plot(        ORBFOL["GCA"]["R"][id1],    ORBFOL["GCA"]["z"][id1],
                colors[4])
    h4.plot(        ORBFOL["GCA"]["R"][id2],    ORBFOL["GCA"]["z"][id2],
                colors[5])

    #**************************************************************************#
    #*                 Finalize and print and show the figure                  #
    #*                                                                         #
    #**************************************************************************#

    h1.set_xlim(0, 5e-6)
    h1.xaxis.set(ticks=[0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6], ticklabels=[])
    h1.yaxis.set(ticks=np.array([-6,-4, -2, 0, 2])*1e-11,
                 ticklabels=[-6, '', '', 0, 2])
    h1.tick_params(axis='y', direction='out')
    h1.tick_params(axis='x', direction='out')
    h1.spines['right'].set_visible(False)
    h1.spines['top'].set_visible(False)
    h1.yaxis.set_ticks_position('left')
    h1.xaxis.set_ticks_position('bottom')
    h1.set(ylabel=r"$\Delta E/E_0\;[10^{-11}]$")

    h2.set_xlim(0, 5e-6)
    h2.xaxis.set(ticks=[0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6], ticklabels=[])
    h2.yaxis.set(ticks=np.array([-4, 0, 4])*1e-3, ticklabels=[-4, 0, 4])
    h2.tick_params(axis='y', direction='out')
    h2.tick_params(axis='x', direction='out')
    h2.spines['right'].set_visible(False)
    h2.spines['top'].set_visible(False)
    h2.yaxis.set_ticks_position('left')
    h2.xaxis.set_ticks_position('bottom')
    h2.set(ylabel=r"$\Delta \mu/\mu_0\;[10^{-3}]$")

    h3.set_xlim(0, 5e-6)
    h3.xaxis.set(ticks=[0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6],
                 ticklabels=[0, 1, 2, 3, 4, 5])
    h3.yaxis.set(ticks=np.array([-6, 0, 6])*1e-4, ticklabels=[-6, 0, 6])
    h3.tick_params(axis='y', direction='out')
    h3.tick_params(axis='x', direction='out')
    h3.spines['right'].set_visible(False)
    h3.spines['top'].set_visible(False)
    h3.yaxis.set_ticks_position('left')
    h3.xaxis.set_ticks_position('bottom')
    h3.set(ylabel=r"$\Delta P/P_0\;[10^{-4}]$", xlabel=r"Time [$10^{-6}$ s]")

    h4.axis('scaled')
    h4.xaxis.set(ticks=[5, 6.5, 8])
    h4.yaxis.set(ticks=[-1.5, 0, 1.5])
    h4.tick_params(axis='y', direction='out')
    h4.tick_params(axis='x', direction='out')
    h4.set(xlabel="$R$ [m]", ylabel="$z$ [m]")

    legend = [r"GO-p", r"GCF-p", r"GCA-p", r"GO-b", r"GCF-b", r"GCA-b"]
    h4.text(5.0, 2.5, legend[0], fontsize=9, color=colors[0])
    h4.text(6.0, 2.5, legend[1], fontsize=9, color=colors[2])
    h4.text(7.0, 2.5, legend[2], fontsize=9, color=colors[4])
    h4.text(5.0, 2.0, legend[3], fontsize=9, color=colors[1])
    h4.text(6.0, 2.0, legend[4], fontsize=9, color=colors[3])
    h4.text(7.0, 2.0, legend[5], fontsize=9, color=colors[5])

    plt.savefig("test_orbitfollowing.png", dpi=300)
    plt.show()

def plot_relerr(axis, x, y, color):
    axis.plot(x, y/y[0] - 1, color)


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
