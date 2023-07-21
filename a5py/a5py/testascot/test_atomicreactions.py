#!/usr/bin/env python3

"""
Test mean free paths for atomic reactions.

This file initializes, runs, and plots test cases for checking that the
ASCOT5 atomic reactions module reproduces analytical mean free paths for
charge-exchange (CX) neutralization and ionization in a plasma.

The tests are done using deuterons and deuterium atoms as markers. We
include orbit following to measure the distance of the path travelled
before a reaction occurs.

Similar tests were reported in
P Ollus et al 2022 Plasma Phys. Control. Fusion 64 035014.

NOTE: This script was originally written with test_coulombcollisions.py as
      a template, i.e., by editing a copy of test_coulombcollisions.py.

File: test_atomicreactions.py
"""


import sys
import copy

import numpy                   as np
import unyt
import scipy

import a5py.ascot5io.options   as options
import a5py.ascot5io.B_GS      as B_GS
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.N0_1D     as N0_1D
import a5py.ascot5io.mrk_gc    as mrk_gc
import a5py.ascot5io.mrk_prt   as mrk_prt
import a5py.ascot5io.boozer    as boozer
import a5py.ascot5io.mhd       as mhd

import a5py.testascot.helpers as helpers
import a5py.physlib.gamma as get_gamma

from a5py.preprocessing.analyticequilibrium import psi0 as psifun
from a5py.preprocessing.asigma_loc_read2hdf5 import read2hdf5 as \
    asigma_read2hdf5

from a5py.ascot5io.ascot5 import Ascot
from a5py.physlib import e, m_e, m_p, c, eps_0


simtime = 1e-1

Te  = 1e3
T0  = 1e3
ne  = 1e20
n0  = 1e16
Eion = 1e5
Entl = 1e5

R0 = 6.2
z0 = 0
Bphi0 = 5.3
psi_mult = 200

# ITER-like but circular equilibrium
psi_coeff = np.array([ 8.629e-02,  3.279e-01,  5.268e-01, -2.366e-01,
                       3.825e-01, -3.573e-01, -1.484e-02,  1.506e-01,
                       7.428e-01, -4.447e-01, -1.084e-01,  1.281e-02, -0.155])

m_Dn = scipy.constants.physical_constants["deuteron mass"][0] * unyt.kg
# An attempt at the exact mass of a deuterium atom
m_Dm = 3.34449439e-27 * unyt.kg


def init():
    """
    Initialize tests

    This function initializes two tests: one that simulates the
    CX neutralization of fast deuterons (TEST_CX), and one that simulates the
    ionization of fast deuterium atoms (TEST_IONZ).

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*              Generate options for TEST_CX and TEST_IONZ                 #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    # Shared simulation settings
    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = simtime
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_COULOMB_COLLISIONS"] = 0
    odict["ENABLE_ATOMIC"]             = 1
    odict["ENABLE_DIST_5D"]            = 1

    # Shared distribution configurations
    odict["DIST_MIN_R"]     = 4
    odict["DIST_MAX_R"]     = 10
    odict["DIST_NBIN_R"]    = 1
    odict["DIST_MIN_PHI"]   = 0
    odict["DIST_MAX_PHI"]   = 360
    odict["DIST_NBIN_PHI"]  = 1
    odict["DIST_MIN_Z"]     = -5
    odict["DIST_MAX_Z"]     = 5
    odict["DIST_NBIN_Z"]    = 1
    odict["DIST_NBIN_PPA"]  = 140
    odict["DIST_NBIN_PPE"]  = 80
    odict["DIST_MIN_TIME"]  = 0
    odict["DIST_MAX_TIME"]  = simtime
    odict["DIST_NBIN_TIME"] = 2

    # Distinct settings for TEST_CX
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-9
    odict["ENDCOND_NEUTRALIZED"]       = 1
    odict["DIST_MIN_PPA"]   = -3.2e6 * m_Dn
    odict["DIST_MAX_PPA"]   =  3.2e6 * m_Dn
    odict["DIST_MIN_PPE"]   = 0      * m_Dn
    odict["DIST_MAX_PPE"]   = 3.2e6  * m_Dn
    options.write_hdf5(helpers.testfn, odict, desc="TEST_CX")

    # Distinct settings for TEST_IONZ
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-13
    odict["ENDCOND_NEUTRALIZED"]       = 0
    odict["ENDCOND_IONIZED"]           = 1
    odict["DIST_MIN_PPA"]   = -3.2e6 * m_Dm
    odict["DIST_MAX_PPA"]   =  3.2e6 * m_Dm
    odict["DIST_MIN_PPE"]   = 0      * m_Dm
    odict["DIST_MAX_PPE"]   = 3.2e6  * m_Dm
    options.write_hdf5(helpers.testfn, odict, desc="TEST_IONZ")

    #**************************************************************************#
    #*       Marker input consisting of deuterons and deuterium atoms          #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 1000
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = (1/Nmrk)* np.ones(ids.shape)
    anum   = 2       * np.ones(ids.shape)
    znum   = 1       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)

    mass   = m_Dn.to("amu") * np.ones(ids.shape)
    pol    = 2 * np.pi      * np.random.rand(1, Nmrk)
    R      = 6.2 + 0.8      * np.cos(pol)
    z      = 0.8            * np.sin(pol)
    charge = 1              * np.ones(ids.shape)
    vtot = get_gamma.vnorm_gamma(get_gamma.gamma_energy(m_Dn,Eion*unyt.eV))
    vr     = -vtot/np.sqrt(2) * np.ones(ids.shape)
    vphi   =  vtot/np.sqrt(2) * np.ones(ids.shape)
    vz     = 0                * np.ones(ids.shape)
    mrk_prt.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, vr, vphi, vz,
                       anum, znum, weight, time, desc="TEST_CX")

    # Place neutral markers closely around the axis with velocity vectors
    # pointing towards the inner-midplane edge to give neutrals enough
    # space to ionize before exiting plasma.
    mass   = m_Dm.to("amu")   * np.ones(ids.shape)
    pol    = 2 * np.pi        * np.random.rand(1, Nmrk)
    R      = 6.2 + 0.1        * np.cos(pol)
    z      = 0.1              * np.sin(pol)
    charge = 0                * np.ones(ids.shape)
    vtot = get_gamma.vnorm_gamma(get_gamma.gamma_energy(m_Dm,Entl*unyt.eV))
    vr     = -vtot/np.sqrt(2) * np.ones(ids.shape)
    vphi   =  vtot/np.sqrt(2) * np.ones(ids.shape)
    vz     = 0                * np.ones(ids.shape)
    mrk_prt.write_hdf5(helpers.testfn, Nmrk, ids, mass,
                       charge, R, phi, z, vr, vphi, vz,
                       anum, znum, weight, time, desc="TEST_IONZ")

    #**************************************************************************#
    #*                 Analytical ITER-like magnetic field                     #
    #*                                                                         #
    #**************************************************************************#
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="TEST_CX")
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="TEST_IONZ")

    #**************************************************************************#
    #*                               Plasma                                    #
    #*                                                                         #
    #**************************************************************************#
    Nrho   = 3
    Nion   = 1
    znum   = np.array([1])
    anum   = np.array([2])
    mass   = np.array([m_Dn.to("amu").v[()]])
    charge = np.array([1])
    rho    = np.array([0, 0.5, 100])
    edens  = ne * np.ones(rho.shape)
    etemp  = Te * np.ones(rho.shape)
    idens  = ne * np.ones((rho.size, Nion))
    itemp  = Te * np.ones(rho.shape)
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, anum, znum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="TEST_CX")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, anum, znum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="TEST_IONZ")

    #**************************************************************************#
    #*                               Neutrals                                  #
    #*                                                                         #
    #**************************************************************************#
    Nntl   = 1
    znum0  = np.array([1])
    anum0  = np.array([2])
    dens0  = n0 * np.ones((rho.size, Nntl))
    temp0  = T0 * np.ones(rho.shape)
    N0_1D.write_hdf5(helpers.testfn, rho[0], rho[-1], Nrho, Nntl, anum0, znum0,
                     dens0, temp0, 1, desc="TEST_CX")

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="TEST_CX")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="TEST_IONZ")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    asigma_dir = "../../../../ascot5/asigma/asigma_loc/"
    for tname in ["TEST_CX", "TEST_IONZ"]:
        W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc=tname)
        boozer.write_hdf5_dummy(helpers.testfn, desc=tname)
        mhd.write_hdf5_dummy(helpers.testfn, desc=tname)
        asigma_read2hdf5(
            fnh5 = helpers.testfn,
            fn_CX_DD0   = asigma_dir + "input.sigmavCX_DD0",
            fn_BMS_H0H  = asigma_dir + "input.BMSsigmav_H0H",
            fn_BMS_H0He = asigma_dir + "input.BMSsigmav_H0He",
            fn_BMS_H0C  = asigma_dir + "input.BMSsigmav_H0C",
            desc = tname)
    N0_1D.write_hdf5_dummy(helpers.testfn, desc="TEST_IONZ")


def run():
    """
    Run tests.
    """
    for test in ["TEST_CX", "TEST_IONZ"]:
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

    a5 = Ascot(helpers.testfn)

    #**************************************************************************#
    #* Read, interpolate and evaluate reaction data for analytical calculation #
    #*                                                                         #
    #**************************************************************************#
    asigma_dir = "../../../../ascot5/asigma/asigma_loc/"
    fn_CX_DD0   = asigma_dir + "input.sigmavCX_DD0"
    fn_BMS_H0H  = asigma_dir + "input.BMSsigmav_H0H"

    # Read CX reaction data of the form sigmav(E, T)
    f = open(fn_CX_DD0, "r")
    if f.mode == 'r':
        lines = f.readlines()
        metadata = lines[0].split(' ')
        N_E = int(metadata[6])
        N_T = int(metadata[7])
        abscissa1 = lines[1].split(' ')
        # Convert the 1st abscissa into floats
        E = np.zeros(N_E)
        for i in range(0, N_E):
            E[i] = float(abscissa1[i])
        abscissa2 = lines[2].split(' ')
        # Convert the 2nd abscissa into floats
        T = np.zeros(N_T)
        for i in range(0, N_T):
            T[i] = float(abscissa2[i])
        ordinate = lines[3].split(' ')
        # Move the rate coefficient data into a 2D array
        sigmavCX_DD0 = np.zeros((N_E,N_T))
        for i_E in range(0, N_E):
            for i_T in range(0, N_T):
                sigmavCX_DD0[i_E,i_T] = float(ordinate[i_T*N_E + i_E])
        # Initialize cubic bivariate spline
        sigmavCX_DD0_interp = scipy.interpolate.RectBivariateSpline(
            E, T, sigmavCX_DD0)

    # Read BMS reaction data of the form BMSsigmav(E, n, T)
    # NOTE: Some of the variables from above are reused
    f = open(fn_BMS_H0H, "r")
    if f.mode == 'r':
        lines = f.readlines()
        metadata = lines[0].split(' ')
        N_E = int(metadata[6])
        N_n = int(metadata[7])
        N_T = int(metadata[8])
        abscissa1 = lines[1].split(' ')
        # Convert the 1st abscissa into floats
        E = np.zeros(N_E)
        for i in range(0, N_E):
            E[i] = float(abscissa1[i])
        abscissa2 = lines[2].split(' ')
        # Convert the 2nd abscissa into floats
        n = np.zeros(N_n)
        for i in range(0, N_n):
            # Convert density unit from cm^-3 to m^-3
            n[i] = 1e6*float(abscissa2[i])
        abscissa3 = lines[3].split(' ')
        # Convert the 3rd abscissa into floats
        T = np.zeros(N_T)
        for i in range(0, N_T):
            T[i] = float(abscissa3[i])
        ordinate = lines[4].split(' ')
        # Move the BMS coefficient data into a 3D array
        BMSsigmav_H0H = np.zeros((N_E, N_n, N_T))
        for i_E in range(0, N_E):
            for i_n in range(0,N_n):
                for i_T in range(0,N_T):
                    # Convert BMS coefficient unit from cm^3s^-1 to m^3s^-1
                    BMSsigmav_H0H[i_E, i_n, i_T] = 1e-6*float(ordinate[
                        i_T*N_n*N_E + i_n*N_E + i_E])
        # Python does not initialize cubic trivariate spline function
        # as it does bivariate spline functions, instead it directly
        # evaluates the target values. Hence, the 3D interpolation is
        # done in a single step below.

    #**************************************************************************#
    #*                      Calculate analytical results                       #
    #*                                                                         #
    #**************************************************************************#
    # The analytical mean free path is calculated as
    #   d_MFP = u/nu = u/(<sigma*v>*n),
    # where u is the fast-particle speed , nu the reaction frequency,
    # <sigma*v> the reaction rate coefficients, and n  the reaction
    # counterpart density.
    u_ion = get_gamma.vnorm_gamma(get_gamma.gamma_energy(m_Dn,Eion*unyt.eV))
    u_ntl = get_gamma.vnorm_gamma(get_gamma.gamma_energy(m_Dm,Entl*unyt.eV))
    sigmav = sigmavCX_DD0_interp(Eion, T0)
    d_MFP_anl_CX   = u_ion/(sigmav*n0)
    # Also calculate analytical mean free time
    t_MFT_anl_CX   = u_ion/u_ion/(sigmav*n0)
    # BMS data energy abscissa in eV/amu. Mass-normalize fast-neutral energy.
    Entl_per_amu = Entl/m_Dm.to("amu")
    BMSsigmav = scipy.interpolate.interpn(
        (E, n, T), BMSsigmav_H0H, np.array([Entl_per_amu, ne, Te]), 'cubic')
    d_MFP_anl_IONZ = u_ntl/(BMSsigmav*ne)
    # Also calculate analytical mean free time
    t_MFT_anl_IONZ = u_ntl/u_ntl/(BMSsigmav*ne)

    #**************************************************************************#
    #*                      Calculate simulated results                        #
    #*                                                                         #
    #**************************************************************************#
    ini_CX   = a5.TEST_CX.inistate.read()
    end_CX   = a5.TEST_CX.endstate.read()
    ini_IONZ = a5.TEST_IONZ.inistate.read()
    end_IONZ = a5.TEST_IONZ.endstate.read()

    # Extract termination times, i.e. time taken before a reaction
    t_CX   = end_CX["time"]
    t_IONZ = end_IONZ["time"]

    # Get number of markers for each test
    Nmrk_CX   = np.size(t_CX)
    Nmrk_IONZ = np.size(t_IONZ)

    # Extract speeds and check that they did not change
    # TEST_CX
    v_CX     = np.zeros(Nmrk_CX)
    v_CX_end = np.zeros(Nmrk_CX)
    for i in range(0,Nmrk_CX):
        # Inistate
        pvec = ( unyt.me*unyt.c/(unyt.me.v[()]*unyt.c.v[()])*
                 np.array([[ini_CX["prprt"][i]],
                           [ini_CX["pphiprt"][i]],
                           [ini_CX["pzprt"][i]]]) )
        vvec = get_gamma.velocity_momentum(ini_CX["mass"][i]*unyt.amu,pvec)
        v_CX[i] = unyt.array.unorm(vvec)
        # Endstate
        pvec = ( unyt.me*unyt.c/(unyt.me.v[()]*unyt.c.v[()])*
                 np.array([[end_CX["prprt"][i]],
                           [end_CX["pphiprt"][i]],
                           [end_CX["pzprt"][i]]]) )
        vvec = get_gamma.velocity_momentum(end_CX["mass"][i]*unyt.amu,pvec)
        v_CX_end[i] = unyt.array.unorm(vvec)
    v_CX     = v_CX*unyt.c/(unyt.c.v[()])
    v_CX_end = v_CX_end*unyt.c/(unyt.c.v[()])
    # TEST_IONZ
    v_IONZ     = np.zeros(Nmrk_IONZ)
    v_IONZ_end = np.zeros(Nmrk_IONZ)
    for i in range(0,Nmrk_IONZ):
        # Inistate
        pvec = ( unyt.me*unyt.c/(unyt.me.v[()]*unyt.c.v[()])*
                 np.array([[ini_IONZ["prprt"][i]],
                           [ini_IONZ["pphiprt"][i]],
                           [ini_IONZ["pzprt"][i]]]) )
        vvec = get_gamma.velocity_momentum(ini_IONZ["mass"][i]*unyt.amu,pvec)
        v_IONZ[i] = unyt.array.unorm(vvec)
        # Endstate
        pvec = ( unyt.me*unyt.c/(unyt.me.v[()]*unyt.c.v[()])*
                 np.array([[end_IONZ["prprt"][i]],
                           [end_IONZ["pphiprt"][i]],
                           [end_IONZ["pzprt"][i]]]) )
        vvec = get_gamma.velocity_momentum(end_IONZ["mass"][i]*unyt.amu,pvec)
        v_IONZ_end[i] = unyt.array.unorm(vvec)
    v_IONZ     = v_IONZ*unyt.c/(unyt.c.v[()])
    v_IONZ_end = v_IONZ_end*unyt.c/(unyt.c.v[()])
    # Check conservation
    if((np.abs(v_CX_end/v_CX - 1) > 1e-6).any()):
        raise(Exception("ERROR, TEST_CX: Marker speed(s) changed. Exiting."))
    if((np.abs(v_IONZ_end/v_IONZ - 1) > 1e-6).any()):
        raise(Exception("ERROR, TEST_IONZ: Marker speed(s) changed. Exiting."))

    # Calculate distances travelled before reaction, i.e. Free Paths, and
    # their mean, i.e. Mean Free Path, for each test
    d_FP_sim_CX    = v_CX*t_CX
    d_FP_sim_IONZ  = v_IONZ*t_IONZ
    d_MFP_sim_CX   = np.mean(d_FP_sim_CX)
    d_MFP_sim_IONZ = np.mean(d_FP_sim_IONZ)
    # Also get Free Times and calculate Mean Free Time (MFT) for each test
    t_FT_sim_CX    = v_CX/v_CX*t_CX
    t_FT_sim_IONZ  = v_IONZ/v_IONZ*t_IONZ
    t_MFT_sim_CX   = np.mean(t_FT_sim_CX)
    t_MFT_sim_IONZ = np.mean(t_FT_sim_IONZ)

    # Estimate standard error of the mean
    err_d_MFP_sim_CX   = np.sqrt(np.sum((d_FP_sim_CX-d_MFP_sim_CX)**2))/Nmrk_CX
    err_d_MFP_sim_IONZ = np.sqrt(
        np.sum((d_FP_sim_IONZ-d_MFP_sim_IONZ)**2))/Nmrk_IONZ
    # Also for the times
    err_t_MFT_sim_CX   = np.sqrt(np.sum((t_FT_sim_CX-t_MFT_sim_CX)**2))/Nmrk_CX
    err_t_MFT_sim_IONZ = np.sqrt(
        np.sum((t_FT_sim_IONZ-t_MFT_sim_IONZ)**2))/Nmrk_IONZ

    #**************************************************************************#
    #*                                Print                                    #
    #*                                                                         #
    #**************************************************************************#
    print(f"Mean free paths (simulated | analytical):\n"
          "  CX neutralization: {:.3f} \u00B1 {:.3f} m | {:.3f} m\n"
          "  ionization:        {:.3f} \u00B1 {:.3f} mm | {:.3f} mm"
          .format(d_MFP_sim_CX.v[()], err_d_MFP_sim_CX.v[()],
                  d_MFP_anl_CX.v[()][0][0],
                  1e3*d_MFP_sim_IONZ.v[()], 1e3*err_d_MFP_sim_IONZ.v[()],
                  1e3*d_MFP_anl_IONZ.v[()][0]))
    print(f"Mean free times (simulated | analytical):\n"
          "  CX neutralization: {:.3f} \u00B1 {:.3f} mus | {:.3f} mus\n"
          "  ionization:        {:.3f} \u00B1 {:.3f} ns | {:.3f} ns\n"
          .format(1e6*t_MFT_sim_CX.v[()], 1e6*err_t_MFT_sim_CX.v[()],
                  1e6*t_MFT_anl_CX.v[()][0][0],
                  1e9*t_MFT_sim_IONZ.v[()], 1e9*err_t_MFT_sim_IONZ.v[()],
                  1e9*t_MFT_anl_IONZ.v[()][0]))


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
