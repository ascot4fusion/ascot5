#!/usr/bin/env python3

"""
Test collisional slowing-down and equilibrium.

This file initializes, runs, and plots test case for checking that ASCOT5
collision operator produces correct distribution for fast-particle slowing-down
and thermal equilibrium.

The tests are done using protons as markers. Even though these tests only
require collision operator, we also include orbit following to be sure that
this does not bias the results.

File: test_coulombcollisions.py
"""

import sys
import copy

import numpy                   as np
import scipy.constants         as constants
import scipy
import matplotlib.pyplot       as plt
import matplotlib.lines        as mlines

import a5py.ascot5io.ascot5    as ascot5
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
m_p_AMU = constants.physical_constants["proton mass in u"][0]
m_p     = constants.physical_constants["proton mass"][0]
m_a_AMU = constants.physical_constants["alpha particle mass in u"][0]
m_a     = constants.physical_constants["alpha particle mass"][0]
m_e_AMU = constants.physical_constants["electron mass in u"][0]
m_e     = constants.physical_constants["electron mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]
eps0    = constants.physical_constants["electric constant"][0]

simtime_th = 2e-2

Te  = 1e3
ne  = 1e20
Eth = 1e3
Esd = 3.5e6

R0 = 6.2
z0 = 0
Bphi0 = 5.3
psi_mult = 200

# ITER-like but circular equilibrium
psi_coeff = np.array([ 8.629e-02,  3.279e-01,  5.268e-01, -2.366e-01,
                       3.825e-01, -3.573e-01, -1.484e-02,  1.506e-01,
                       7.428e-01, -4.447e-01, -1.084e-01,  1.281e-02, -0.155])

def init():
    """
    Initialize tests

    This function initializes three tests that simulate a thermal equilibrium
    in different modes: THERMAL_GO, THERMAL_GCF, and THERMAL_GCA. Three
    additional simulations are initialized that simulate alpha particle
    slowing-down in different modes: SLOWING_GO, SLOWING_GCF, and SLOWING_GCA.

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*       Generate options for THERMAL_GO, THERMAL_GCF, and THERMAL_GCA     #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = simtime_th
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["ENABLE_DIST_5D"]            = 1

    odict["DIST_MIN_R"]     = 4
    odict["DIST_MAX_R"]     = 10
    odict["DIST_NBIN_R"]    = 1
    odict["DIST_MIN_PHI"]   = 0
    odict["DIST_MAX_PHI"]   = 360
    odict["DIST_NBIN_PHI"]  = 1
    odict["DIST_MIN_Z"]     = -5
    odict["DIST_MAX_Z"]     = 5
    odict["DIST_NBIN_Z"]    = 1
    odict["DIST_MIN_VPA"]   = -1.5e6
    odict["DIST_MAX_VPA"]   =  1.5e6
    odict["DIST_NBIN_VPA"]  = 140
    odict["DIST_MIN_VPE"]   = 0
    odict["DIST_MAX_VPE"]   = 1.5e6
    odict["DIST_NBIN_VPE"]  = 80
    odict["DIST_MIN_TIME"]  = 0
    odict["DIST_MAX_TIME"]  = simtime_th
    odict["DIST_NBIN_TIME"] = 2

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    options.write_hdf5(helpers.testfn, odict, desc="THERMAL_GO")

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 2e-8
    options.write_hdf5(helpers.testfn, odict, desc="THERMAL_GCF")

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 1e-6
    odict["ADAPTIVE_TOL_COL"]          = 1e-2
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    options.write_hdf5(helpers.testfn, odict, desc="THERMAL_GCA")


    #**************************************************************************#
    #*  Generate options for SLOWING_GO, SLOWING_GCF, and SLOWING_GCA          #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["ENDCOND_ENERGYLIM"]         = 1
    odict["ENDCOND_MIN_ENERGY"]        = 50e3
    odict["ENDCOND_MIN_THERMAL"]       = 0
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["ENABLE_DIST_5D"]            = 1

    odict["DIST_MIN_R"]    = 4
    odict["DIST_MAX_R"]    = 10
    odict["DIST_NBIN_R"]   = 1
    odict["DIST_MIN_PHI"]  = 0
    odict["DIST_MAX_PHI"]  = 360
    odict["DIST_NBIN_PHI"] = 1
    odict["DIST_MIN_Z"]    = -5
    odict["DIST_MAX_Z"]    = 5
    odict["DIST_NBIN_Z"]   = 1
    odict["DIST_MIN_VPA"]  = -2e7
    odict["DIST_MAX_VPA"]  =  2e7
    odict["DIST_NBIN_VPA"] = 200
    odict["DIST_MIN_VPE"]  = 0
    odict["DIST_MAX_VPE"]  = 2e7
    odict["DIST_NBIN_VPE"] = 100

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 2e-9
    options.write_hdf5(helpers.testfn, odict, desc="SLOWING_GO")

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 3e-8
    options.write_hdf5(helpers.testfn, odict, desc="SLOWING_GCF")

    odict["SIM_MODE"]                  = 2
    odict["ENABLE_ADAPTIVE"]           = 1
    odict["ADAPTIVE_TOL_ORBIT"]        = 1e-6
    odict["ADAPTIVE_TOL_COL"]          = 1e-2
    odict["ADAPTIVE_MAX_DRHO"]         = 0.1
    odict["ADAPTIVE_MAX_DPHI"]         = 10
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
    options.write_hdf5(helpers.testfn, odict, desc="SLOWING_GCA")

    #**************************************************************************#
    #*            Marker input consisting of protons and alfas                 #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 20
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = (1/Nmrk)* np.ones(ids.shape)
    mass   = m_p_AMU * np.ones(ids.shape)
    charge = 1       * np.ones(ids.shape)
    anum   = 1       * np.ones(ids.shape)
    znum   = 1       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    zeta   = 2 * np.pi * np.random.rand(1,Nmrk)
    energy = Eth     * np.ones(ids.shape)
    pitch  = 0.5     * np.ones(ids.shape)

    pol    = 2*np.pi*np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8 * np.cos(pol)
    z      = 0.8 * np.sin(pol)
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="THERMAL_GO")

    pol    = 2*np.pi*np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8 * np.cos(pol)
    z      = 0.8 * np.sin(pol)
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="THERMAL_GCF")

    pol    = 2*np.pi*np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8 * np.cos(pol)
    z      = 0.8 * np.sin(pol)
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="THERMAL_GCA")

    Nmrk   = 200
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = (1/Nmrk)* np.ones(ids.shape)
    mass   = m_a_AMU * np.ones(ids.shape)
    charge = 2       * np.ones(ids.shape)
    anum   = 4       * np.ones(ids.shape)
    znum   = 2       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    zeta   = 2 * np.pi * np.random.rand(1,Nmrk)
    energy = Esd * np.ones(ids.shape)

    pol    = 2*np.pi*np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8 * np.cos(pol)
    z      = 0.8 * np.sin(pol)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="SLOWING_GO")

    pol    = 2*np.pi*np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8 * np.cos(pol)
    z      = 0.8 * np.sin(pol)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="SLOWING_GCF")

    pol    = 2*np.pi*np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8 * np.cos(pol)
    z      = 0.8 * np.sin(pol)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="SLOWING_GCA")

    #**************************************************************************#
    #*                 Analytical ITER-like magnetic field                     #
    #*                                                                         #
    #**************************************************************************#
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="THERMAL_GO")
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="THERMAL_GCF")
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="THERMAL_GCA")
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="SLOWING_GO")
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="SLOWING_GCF")
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="SLOWING_GCA")

    #**************************************************************************#
    #*                               Plasma                                    #
    #*                                                                         #
    #**************************************************************************#
    Nrho   = 3
    Nion   = 1
    znum   = np.array([1])
    anum   = np.array([1])
    mass   = np.array([1])
    charge = np.array([1])
    rho    = np.array([0, 0.5, 100])
    edens  = ne  * np.ones(rho.shape)
    etemp  = Te  * np.ones(rho.shape)
    idens  = ne  * np.ones((rho.size, Nion))
    itemp  = 1e3 * np.ones(rho.shape)
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="THERMAL_GO")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="THERMAL_GCF")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="THERMAL_GCA")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="SLOWING_GO")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="SLOWING_GCF")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="SLOWING_GCA")

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="THERMAL_GO")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="THERMAL_GCF")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="THERMAL_GCA")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="SLOWING_GO")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="SLOWING_GCF")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="SLOWING_GCA")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="THERMAL_GO")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="THERMAL_GCF")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="THERMAL_GCA")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="SLOWING_GO")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="SLOWING_GCF")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="SLOWING_GCA")

    N0_3D.write_hdf5_dummy(helpers.testfn, desc="THERMAL_GO")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="THERMAL_GCF")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="THERMAL_GCA")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="SLOWING_GO")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="SLOWING_GCF")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="SLOWING_GCA")

def run():
    """
    Run tests.
    """
    for test in ["THERMAL_GO", "THERMAL_GCF", "THERMAL_GCA", \
                 "SLOWING_GO", "SLOWING_GCF", "SLOWING_GCA"]:
        helpers.set_and_run(test)

def check():
    """
    Plot the results of these tests.

    This function makes four plots.
    - One that shows marker equilibrium energy distributions (for all modes) and
      compares them with a Maxwellian distribution.
    - One that shows the pitch distributions in equilibrium.
    - One that shows marker slowing-down energy distributions and compares them
      to the analytical result.
    - One that shows pitch distribution during the slowing-down.
    """

    a5 = ascot5.Ascot(helpers.testfn)

    ts_GO  = np.mean(a5.SLOWING_GO.endstate["time"])
    ts_GCF = np.mean(a5.SLOWING_GCF.endstate["time"])
    ts_GCA = np.mean(a5.SLOWING_GCA.endstate["time"])

    #**************************************************************************#
    #*      Initialize histogram grids and evaluate analytical results         #
    #*                                                                         #
    #**************************************************************************#
    THERMAL = {}
    THERMAL["Egrid"] = np.linspace( 0,   4e4, 200)
    THERMAL["xigrid"]  = np.linspace(-1,     1, 50)
    SLOWING = {}
    SLOWING["Egrid"] = np.linspace( 1, 4.0e6, 250)
    SLOWING["xigrid"]  = np.linspace(-1,     1, 20)

    alphaZ = 2
    clog   = 16
    vth    = np.sqrt(2*Te*e / m_e)
    vcrit  = vth * np.power( (3.0*np.sqrt(np.pi)/4.0) * (m_e / m_p) , 1/3.0)
    Ecrit  = 0.5 * m_a * vcrit * vcrit / e
    ts     = 3 * np.sqrt( np.power(2*np.pi * Te * e, 3) / m_e ) * eps0 * eps0 \
             * m_a /( alphaZ * alphaZ * np.power(e, 4) * ne * clog)

    heaviside = np.logical_and(SLOWING["Egrid"] <= Esd,
                               SLOWING["Egrid"] >= 50*Te)

    THERMAL["analytical"]  = 2 * np.sqrt(THERMAL["Egrid"]/np.pi) \
                             * np.power(Te,-3.0/2) \
                             * np.exp(-THERMAL["Egrid"]/Te)
    THERMAL["analytical"] *= (simtime_th/2)
    SLOWING["analytical"] = heaviside * ts \
                            / ( ( 1 + np.power(Ecrit/SLOWING["Egrid"], 3.0/2) )\
                                * 2 * SLOWING["Egrid"] )

    # ts is slowing down rate which gives the slowing down time as
    # t_sd = ts*log(v_0 / v_th) = 0.5*ts*log(E_0/E_th)
    slowingdowntime = 0.5*ts*np.log(Esd/(50*Te))

    #**************************************************************************#
    #*            Evaluate thermal distributions in energy and pitch           #
    #*                                                                         #
    #**************************************************************************#

    for mode in ["GO", "GCF", "GCA"]:
        distobj = a5["THERMAL_" + mode]["dist5d"]
        Exidist = distobj.get_E_xi_dist(E_edges=THERMAL["Egrid"],
                                        xi_edges=THERMAL["xigrid"],
                                        r=0, phi=0, z=0, charge=0, time=[1, 2])

        Edist  = Exidist
        xidist = copy.deepcopy(Exidist)
        distobj.get_dist(dist=Edist, pitch=0)
        distobj.get_dist(dist=xidist, energy=0)

        THERMAL[mode + "_Edist"]  = Edist["distribution"]
        THERMAL[mode + "_xidist"] = xidist["distribution"]
        THERMAL["E"]  = Edist["energy"]
        THERMAL["xi"] = xidist["pitch"]

    #**************************************************************************#
    #*        Evaluate slowing-down distributions in energy and pitch          #
    #*                                                                         #
    #**************************************************************************#

    for mode in ["GO", "GCF", "GCA"]:
        distobj = a5["SLOWING_" + mode]["dist5d"]

        Exidist = distobj.get_E_xi_dist(E_edges=SLOWING["Egrid"],
                                        xi_edges=SLOWING["xigrid"],
                                        r=0, phi=0, z=0, charge=0, time=0)

        Edist  = Exidist
        xidist = copy.deepcopy(Exidist)
        distobj.get_dist(dist=Edist,   pitch=0)
        distobj.get_dist(dist=xidist, energy=0)

        SLOWING[mode + "_Edist"]  = Edist["distribution"]
        SLOWING[mode + "_xidist"] = xidist["distribution"]
        SLOWING["E"]  = Edist["energy"]
        SLOWING["xi"] = xidist["pitch"]

    #**************************************************************************#
    #*                                 Plot                                    #
    #*                                                                         #
    #**************************************************************************#

    cy = plt.rcParams['axes.prop_cycle'].by_key()['color']

    f = plt.figure(figsize=(11.9/2.54, 8/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    c = ['b', 'g', 'r']
    h1 = f.add_subplot(1,4,1)
    h1.set_position([0.12, 0.67, 0.45, 0.3], which='both')

    h2 = f.add_subplot(1,4,2)
    h2.set_position([0.72, 0.67, 0.25, 0.3], which='both')

    h3 = f.add_subplot(1,4,3)
    h3.set_position([0.12, 0.17, 0.45, 0.3], which='both')

    h4 = f.add_subplot(1,4,4)
    h4.set_position([0.72, 0.17, 0.25, 0.3], which='both')

    h1.plot(THERMAL["Egrid"], THERMAL["analytical"], 'black')
    h1.plot(THERMAL["E"], THERMAL["GO_Edist"],
            marker='*', markersize=8, markevery=3, linestyle='none', alpha=0.5)
    h1.plot(THERMAL["E"], THERMAL["GCF_Edist"],
            marker='.', markersize=10, markevery=3, linestyle='none', alpha=0.5)
    h1.plot(THERMAL["E"], THERMAL["GCA_Edist"],
            marker='^', markersize=5, markevery=3, linestyle='none', alpha=0.5)

    h2.plot(np.array([-1, 1]), np.array([0.5, 0.5])*simtime_th/2, 'black')
    h2.plot(THERMAL["xi"], THERMAL["GO_xidist"], alpha=0.5)
    h2.plot(THERMAL["xi"], THERMAL["GCF_xidist"], alpha=0.5)
    h2.plot(THERMAL["xi"], THERMAL["GCA_xidist"], alpha=0.5)

    dense = 16
    h3.plot(SLOWING["Egrid"], SLOWING["analytical"], 'black')
    h3.plot(SLOWING["E"][:dense], SLOWING["GO_Edist"][:dense], marker='*',
            markersize=8, markevery=2, linestyle='none', alpha=0.5, color=cy[0])
    h3.plot(SLOWING["E"][dense:], SLOWING["GO_Edist"][dense:], marker='*',
            markersize=8, markevery=20, linestyle='none', alpha=0.5,color=cy[0])
    h3.plot(SLOWING["E"][:dense], SLOWING["GCF_Edist"][:dense], marker='.',
            markersize=10, markevery=2, linestyle='none', alpha=0.5,color=cy[1])
    h3.plot(SLOWING["E"][dense:], SLOWING["GCF_Edist"][dense:], marker='.',
            markersize=10, markevery=20, linestyle='none', alpha=0.5,color=cy[1])
    h3.plot(SLOWING["E"][:dense], SLOWING["GCA_Edist"][:dense], marker='^',
            markersize=5, markevery=2, linestyle='none', alpha=0.5, color=cy[2])
    h3.plot(SLOWING["E"][dense:], SLOWING["GCA_Edist"][dense:], marker='^',
            markersize=5, markevery=20, linestyle='none', alpha=0.5, color=cy[2])

    h4.plot(np.array([-1, 1]), np.array([0.5, 0.5])*slowingdowntime, 'black')
    h4.plot(SLOWING["xi"], SLOWING["GO_xidist"], alpha=0.5)
    h4.plot(SLOWING["xi"], SLOWING["GCF_xidist"], alpha=0.5)
    h4.plot(SLOWING["xi"], SLOWING["GCA_xidist"], alpha=0.5)

    #**************************************************************************#
    #*                         Decorate and finish                             #
    #*                                                                         #
    #**************************************************************************#
    h1.set_xlim(0, 1e4)
    h1.xaxis.set(ticks=[0, 2.5e3, 5e3, 7.5e3, 1e4],
                 ticklabels=['0', '2.5', '5', '7.5', '10'])
    h1.yaxis.set(ticks=np.array([0, 2, 4, 6])*1e-6,
                 ticklabels=['0', '2', '4', '6'])
    h1.tick_params(axis='y', direction='out')
    h1.tick_params(axis='x', direction='out')
    h1.spines['right'].set_visible(False)
    h1.spines['top'].set_visible(False)
    h1.yaxis.set_ticks_position('left')
    h1.xaxis.set_ticks_position('bottom')
    h1.set(ylabel=r"$f(E) \times 10^{-6}$ [1/eV]")
    h1.set(xlabel=r"Energy $\times 10^{3}$ [eV]")

    h2.set_xlim(-1, 1)
    h2.set_ylim(0, 1*simtime_th/2)
    h2.xaxis.set(ticks=[-1, -0.5, 0, 0.5, 1])
    h2.yaxis.set(ticks=np.array([0, 0.5, 1])*simtime_th/2,
                 ticklabels=['0', '5', '10'])
    h2.tick_params(axis='y', direction='out')
    h2.tick_params(axis='x', direction='out')
    h2.spines['right'].set_visible(False)
    h2.spines['top'].set_visible(False)
    h2.yaxis.set_ticks_position('left')
    h2.xaxis.set_ticks_position('bottom')
    h2.set(ylabel=r"$f(\xi) \times 10^{-4}$")
    h2.set(xlabel=r"Pitch")

    h3.set_xlim(0, 4e6)
    h3.xaxis.set(ticks=[0, 1e6, 2e6, 3e6, 4e6],
                 ticklabels=['0', '1', '2', '3', '4'])
    h3.yaxis.set(ticks=np.array([0, 2, 4, 6, 8])*1e-8,
                 ticklabels=['0', '2', '4', '6', '8'])
    h3.tick_params(axis='y', direction='out')
    h3.tick_params(axis='x', direction='out')
    h3.spines['right'].set_visible(False)
    h3.spines['top'].set_visible(False)
    h3.yaxis.set_ticks_position('left')
    h3.xaxis.set_ticks_position('bottom')
    h3.set(ylabel=r"$f(E) \times 10^{-8}$ [1/eV]")
    h3.set(xlabel=r"Energy $\times 10^{6}$ [eV]")

    h4.set_xlim(-1, 1)
    h4.set_ylim(0, 1 * ts)
    h4.xaxis.set(ticks=[-1, -0.5, 0, 0.5, 1])
    h4.yaxis.set(ticks=[0, 0.015, 0.03],
                 ticklabels=['0', '15', '30'])
    h4.tick_params(axis='y', direction='out')
    h4.tick_params(axis='x', direction='out')
    h4.spines['right'].set_visible(False)
    h4.spines['top'].set_visible(False)
    h4.yaxis.set_ticks_position('left')
    h4.xaxis.set_ticks_position('bottom')
    h4.set(ylabel=r"$f(\xi)\times 10^{-3}$")
    h4.set(xlabel=r"Pitch")

    l1 = mlines.Line2D([], [], color='black', linestyle='-',
                       label='Analytical',axes=h1)
    l2 = mlines.Line2D([], [], color=cy[0], marker='*', linestyle='-',
                       markersize=6, label='GO',axes=h1)
    l3 = mlines.Line2D([], [], color=cy[1], marker='.', linestyle='-',
                       markersize=8, label='GCF',axes=h1)
    l4 = mlines.Line2D([], [], color=cy[2], marker='^', linestyle='-',
                       markersize=4, label='GCA',axes=h1)

    h1.legend(handles=[l1, l2, l3, l4], loc='upper right', frameon=False,
              fontsize=8)

    frm  = lambda x: "%3.0f ms" % (x*1e3)
    text = r"$\tau_s$:      " + frm(slowingdowntime) + "\n" \
           + "GO:    " + frm(ts_GO) + "\n" \
           + "GCF:  " + frm(ts_GCF) + "\n" \
           + "GCA: " + frm(ts_GCA)
    h3.text(2e6, 3e-8, text, fontsize=9)

    plt.savefig("test_coulombcollisions.png", dpi=300)
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
