"""
Test collisional slowing-down and equilibrium.

File: test_coulombcollisions.py
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
m_p_AMU = constants.physical_constants["proton mass in u"][0]
m_p     = constants.physical_constants["proton mass"][0]
m_e_AMU = constants.physical_constants["electron mass in u"][0]
m_e     = constants.physical_constants["electron mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

Te   = 1e3
ne   = 1e20

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
    #*                     Generate options for CCOLL-GO                       #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-6
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-3
    odict["ENABLE_ORBIT_FOLLOWING"]    = 0
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"] = 1
    odict["DIST_MIN_R"]    = 4
    odict["DIST_MAX_R"]    = 8
    odict["DIST_NBIN_R"]   = 1
    odict["DIST_MIN_phi"]  = 0
    odict["DIST_MAX_phi"]  = 360
    odict["DIST_NBIN_phi"] = 1
    odict["DIST_MIN_z"]    = 5
    odict["DIST_MAX_z"]    = 5
    odict["DIST_NBIN_z"]   = 1
    odict["DIST_MIN_vpa"]  = -2e6
    odict["DIST_MAX_vpa"]  =  2e6
    odict["DIST_NBIN_vpa"] = 100
    odict["DIST_MIN_vpe"]  = 0
    odict["DIST_MAX_vpe"]  = 2e6
    odict["DIST_NBIN_vpe"] = 50

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="THERMAL_GO")

    odict["ENDCOND_SIMTIMELIM"]               = 0
    odict["ENDCOND_ENERGYLIM"]                = 1
    odict["ENDCOND_MIN_ENERGY"]               = 1e3
    odict["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = 0
    odict["DIST_MIN_vpa"]  = -5e6
    odict["DIST_MAX_vpa"]  =  5e6
    odict["DIST_NBIN_vpa"] = 200
    odict["DIST_MIN_vpe"]  = 0
    odict["DIST_MAX_vpe"]  = 5e6
    odict["DIST_NBIN_vpe"] = 100

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="SLOWING_GO")

    #**************************************************************************#
    #*                     Generate options for CCOLL-GCF                      #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-9
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 0
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"] = 1
    odict["DIST_MIN_R"]    = 4
    odict["DIST_MAX_R"]    = 8
    odict["DIST_NBIN_R"]   = 1
    odict["DIST_MIN_phi"]  = 0
    odict["DIST_MAX_phi"]  = 360
    odict["DIST_NBIN_phi"] = 1
    odict["DIST_MIN_z"]    = 5
    odict["DIST_MAX_z"]    = 5
    odict["DIST_NBIN_z"]   = 1
    odict["DIST_MIN_vpa"]  = -2e6
    odict["DIST_MAX_vpa"]  =  2e6
    odict["DIST_NBIN_vpa"] = 100
    odict["DIST_MIN_vpe"]  = 0
    odict["DIST_MAX_vpe"]  = 2e6
    odict["DIST_NBIN_vpe"] = 50

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="THERMAL_GCF")

    odict["ENDCOND_SIMTIMELIM"]               = 0
    odict["ENDCOND_ENERGYLIM"]                = 1
    odict["ENDCOND_MIN_ENERGY"]               = 1e3
    odict["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = 0
    odict["DIST_MIN_vpa"]  = -5e6
    odict["DIST_MAX_vpa"]  =  5e6
    odict["DIST_NBIN_vpa"] = 200
    odict["DIST_MIN_vpe"]  = 0
    odict["DIST_MAX_vpe"]  = 5e6
    odict["DIST_NBIN_vpe"] = 100

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="SLOWING_GCF")

    #**************************************************************************#
    #*                     Generate options for CCOLL-GCA                      #
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
    odict["ENABLE_ORBIT_FOLLOWING"]    = 0
    odict["ENABLE_COULOMB_COLLISIONS"] = 1
    odict["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"] = 1
    odict["DIST_MIN_R"]    = 4
    odict["DIST_MAX_R"]    = 8
    odict["DIST_NBIN_R"]   = 1
    odict["DIST_MIN_phi"]  = 0
    odict["DIST_MAX_phi"]  = 360
    odict["DIST_NBIN_phi"] = 1
    odict["DIST_MIN_z"]    = 5
    odict["DIST_MAX_z"]    = 5
    odict["DIST_NBIN_z"]   = 1
    odict["DIST_MIN_vpa"]  = -2e6
    odict["DIST_MAX_vpa"]  =  2e6
    odict["DIST_NBIN_vpa"] = 100
    odict["DIST_MIN_vpe"]  = 0
    odict["DIST_MAX_vpe"]  = 2e6
    odict["DIST_NBIN_vpe"] = 50

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="THERMAL_GCA")

    odict["ENDCOND_SIMTIMELIM"]               = 0
    odict["ENDCOND_ENERGYLIM"]                = 1
    odict["ENDCOND_MIN_ENERGY"]               = 1e3
    odict["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = 0
    odict["DIST_MIN_vpa"]  = -5e6
    odict["DIST_MAX_vpa"]  =  5e6
    odict["DIST_NBIN_vpa"] = 200
    odict["DIST_MIN_vpe"]  = 0
    odict["DIST_MAX_vpe"]  = 5e6
    odict["DIST_NBIN_vpe"] = 100

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="SLOWING_GCA")

    #**************************************************************************#
    #*                    Marker input consisting of protons                   #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 20
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = np.ones(ids.shape)
    mass   = m_p_AMU * np.ones(ids.shape)
    charge = 1       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    R      = 8       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    z      = 0       * np.ones(ids.shape)
    theta  = 2 * np.pi * np.random.rand(1,Nmrk)
    energy = 1e3     * np.ones(ids.shape)
    pitch  = 0.9     * np.ones(ids.shape)

    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, theta,
                   weight, time, desc="THERMAL_GO")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, theta,
                   weight, time, desc="THERMAL_GCF")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, theta,
                   weight, time, desc="THERMAL_GCA")

    Nmrk   = 2
    ids    = np.linspace(1,Nmrk,Nmrk)
    weight = np.ones(ids.shape)
    mass   = m_p_AMU * np.ones(ids.shape)
    charge = 1       * np.ones(ids.shape)
    time   = 0       * np.ones(ids.shape)
    R      = 8       * np.ones(ids.shape)
    phi    = 90      * np.ones(ids.shape)
    z      = 0       * np.ones(ids.shape)
    theta  = 2 * np.pi * np.random.rand(1,Nmrk)
    energy = 3.5e6 * np.ones(ids.shape)
    pitch  = 1 - 2 * np.random.rand(1,Nmrk)
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, theta,
                   weight, time, desc="SLOWING_GO")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, theta,
                   weight, time, desc="SLOWING_GCF")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass,
                   charge, R, phi, z, energy, pitch, theta,
                   weight, time, desc="SLOWING_GCA")

    #**************************************************************************#
    #*     Uniform magnetic field with values scanned from Bmin to Bmax        #
    #*                                                                         #
    #**************************************************************************#
    B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="THERMAL_GO")
    B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="THERMAL_GCF")
    B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="THERMAL_GCA")
    B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="SLOWING_GO")
    B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="SLOWING_GCF")
    B_GS.write_hdf5(test_ascot.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="SLOWING_GCA")

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
    idens = ne  * np.ones((rho.size, Nion))
    itemp = 1e3 * np.ones(rho.shape)
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="THERMAL_GO")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="THERMAL_GCF")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="THERMAL_GCA")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="SLOWING_GO")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="SLOWING_GCF")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="SLOWING_GCA")

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="THERMAL_GO")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="THERMAL_GCF")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="THERMAL_GCA")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="SLOWING_GO")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="SLOWING_GCF")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="SLOWING_GCA")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="THERMAL_GO")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="THERMAL_GCF")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="THERMAL_GCA")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="SLOWING_GO")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="SLOWING_GCF")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="SLOWING_GCA")

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
                     desc="THERMAL_GO")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens,
                     desc="THERMAL_GCF")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens,
                     desc="THERMAL_GCA")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens,
                     desc="SLOWING_GO")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens,
                     desc="SLOWING_GCF")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens,
                     desc="SLOWING_GCA")

def run():
    """
    Run tests.
    """
    for test in ["THERMAL_GO", "THERMAL_GCF", "THERMAL_GCA", \
                 "SLOWING_GO", "SLOWING_GCF", "SLOWING_GCA"]:
        test_ascot.set_and_run(test)

def check():

    a5 = ascot5.Ascot(test_ascot.testfn)

    THERMAL = {}
    THERMAL["energygrid"] = np.linspace( 0,   6e3, 50)
    THERMAL["pitchgrid"]  = np.linspace(-1,     1, 50)
    SLOWING = {}
    SLOWING["energygrid"] = np.linspace( 0, 1.1e5, 50)
    SLOWING["pitchgrid"]  = np.linspace(-1,     1, 50)
    SLOWING["timegrid"]   = np.linspace( 0,  0.03, 50)

    THERMAL["GO_energydist"]  = np.histogram(a5["THERMAL_GO"].endstate["Ekin"]/e,
                                             bins=THERMAL["energygrid"],
                                             density=True)[0]
    THERMAL["GO_pitchdist"]   = np.histogram(a5["THERMAL_GO"].endstate["pitch"],
                                             bins=THERMAL["pitchgrid"],
                                             density=True)[0]
    THERMAL["GCF_energydist"] = np.histogram(a5["THERMAL_GCF"].endstate["Ekin"]/e,
                                             bins=THERMAL["energygrid"],
                                             density=True)[0]
    THERMAL["GCF_pitchdist"]  = np.histogram(a5["THERMAL_GCF"].endstate["pitch"],
                                             bins=THERMAL["pitchgrid"],
                                             density=True)[0]
    THERMAL["GCA_energydist"] = np.histogram(a5["THERMAL_GCA"].endstate["Ekin"]/e,
                                             bins=THERMAL["energygrid"],
                                             density=True)[0]
    THERMAL["GCA_pitchdist"]  = np.histogram(a5["THERMAL_GCA"].endstate["pitch"],
                                             bins=THERMAL["pitchgrid"],
                                             density=True)[0]

    SLOWING["GO_pitchdist"] = np.histogram(a5["SLOWING_GO"].endstate["pitch"],
                                           bins=SLOWING["pitchgrid"],
                                           density=True)[0]
    SLOWING["GO_timedist"]  = np.histogram(a5["SLOWING_GO"].endstate["time"],
                                           bins=SLOWING["timegrid"],
                                           density=True)[0]
    SLOWING["GCF_pitchdist"] = np.histogram(a5["SLOWING_GCF"].endstate["pitch"],
                                            bins=SLOWING["pitchgrid"],
                                            density=True)[0]
    SLOWING["GCF_timedist"]  = np.histogram(a5["SLOWING_GCF"].endstate["time"],
                                            bins=SLOWING["timegrid"],
                                            density=True)[0]
    SLOWING["GCA_pitchdist"] = np.histogram(a5["SLOWING_GCA"].endstate["pitch"],
                                            bins=SLOWING["pitchgrid"],
                                            density=True)[0]
    SLOWING["GCA_timedist"]  = np.histogram(a5["SLOWING_GCA"].endstate["time"],
                                            bins=SLOWING["timegrid"],
                                            density=True)[0]

    SLOWING["GO_energypitchdist"] \
        = a5["SLOWING_GO"]["R_phi_z_vpa_vpe_t_q"].get_energypitchdist(SLOWING["energygrid"],
                                               SLOWING["pitchgrid"],
                                               m_p)["ordinate"]
    SLOWING["GCF_energypitchdist"] \
        = a5["SLOWING_GCF"]["R_phi_z_vpa_vpe_t_q"].get_energypitchdist(SLOWING["energygrid"],
                                                SLOWING["pitchgrid"],
                                                m_p)["ordinate"]
    SLOWING["GCA_energypitchdist"] \
        = a5["SLOWING_GCA"]["R_phi_z_vpa_vpe_t_q"].get_energypitchdist(SLOWING["energygrid"],
                                                SLOWING["pitchgrid"],
                                                m_p)["ordinate"]

    print(SLOWING["GO_energypitchdist"].shape)
    SLOWING["GO_energydist"] = np.squeeze(SLOWING["GO_energypitchdist"]).sum(axis=(1))
    print(SLOWING["GO_energydist"].shape)
    #SLOWING["GO_energydist"] = np.sum(  SLOWING["GO_energydist"]
    #                                    * (SLOWING["energygrid"][1]
    #                                       - SLOWING["energygrid"][0]) )

    THERMAL["energygrid"] = ( THERMAL["energygrid"][:-1]
                              + THERMAL["energygrid"][1:] ) / 2
    THERMAL["pitchgrid"]  = ( THERMAL["pitchgrid"][:-1]
                              + THERMAL["pitchgrid"][1:] ) / 2
    SLOWING["energygrid"] = ( SLOWING["energygrid"][:-1]
                              + SLOWING["energygrid"][1:] ) / 2
    SLOWING["pitchgrid"]  = ( SLOWING["pitchgrid"][:-1]
                              + SLOWING["pitchgrid"][1:] ) / 2
    SLOWING["timegrid"]   = ( SLOWING["timegrid"][:-1]
                              + SLOWING["timegrid"][1:] ) / 2

    #thEgrid    = (thEgrid[:-1] + thEgrid[1:]) / 2
    #pitchgrid  = (pitchgrid[:-1] + pitchgrid[1:]) / 2
    #sdEgrid    = (sdEgrid[:-1] + sdEgrid[1:]) / 2
    #sdtimegrid = (sdtimegrid[:-1] + sdtimegrid[1:]) / 2

    #anthEkin   = np.sqrt(thEgrid)*np.exp(-thEgrid/Temp)
    #anthEkin   = anthEkin/np.sum(anthEkin*(thEgrid[1]-thEgrid[0]))

    #anthpitch  = np.ones(pitchgrid.shape)
    #anthpitch  = anthpitch/np.sum(anthpitch*(pitchgrid[1]-pitchgrid[0]))
    #ansdpitch  = anthpitch

    #v0         = E2v(3.5e6*e, m_p)
    #v          = E2v(sdEgrid*e, m_p)
    #vth        = E2v(Temp*e, m_e)
    #vcrit      = vth * np.power( (3.0*np.sqrt(np.pi)/4.0) * m_e / m_p , 1/3.0)
    #ansdEkin   = heaviside(v0-v, 1) *  v / (np.power(vcrit,3) + np.power(v,3))
    #ansdEkin   = ansdEkin/np.sum(ansdEkin*(sdEgrid[1]-sdEgrid[0]))




    c = ['b', 'g', 'r']
    f  = plt.figure()
    a1 = f.add_subplot(3,3,1)
    a2 = f.add_subplot(3,3,2)
    a3 = f.add_subplot(3,3,3)
    a4 = f.add_subplot(3,3,4)
    a5 = f.add_subplot(3,3,5)
    a6 = f.add_subplot(3,3,6)
    a7 = f.add_subplot(3,3,7)
    a8 = f.add_subplot(3,3,8)
    a9 = f.add_subplot(3,3,9)

    #a1.plot(thEgrid,anthEkin)
    #a2.plot(pitchgrid,anthpitch)

    #a3.plot(sdEgrid,ansdEkin)
    a4.plot(0,0)
    #a5.plot(pitchgrid,ansdpitch)

    a1.plot(THERMAL["energygrid"], THERMAL["GO_energydist"])
    a1.plot(THERMAL["energygrid"], THERMAL["GCF_energydist"])
    a1.plot(THERMAL["energygrid"], THERMAL["GCA_energydist"])

    a2.plot(THERMAL["pitchgrid"], THERMAL["GO_pitchdist"])
    a2.plot(THERMAL["pitchgrid"], THERMAL["GCF_pitchdist"])
    a2.plot(THERMAL["pitchgrid"], THERMAL["GCA_pitchdist"])

    print(SLOWING["GO_energydist"])
    a3.plot(SLOWING["energygrid"], SLOWING["GO_energydist"])
    a4.plot(SLOWING["timegrid"], SLOWING["GO_timedist"])
    a5.plot(SLOWING["pitchgrid"], SLOWING["GO_pitchdist"])
    plt.show()
    return

    for i in range(0,3):
        a3.plot(sdEgrid,sdEdist[i])
        a4.plot(sdtimegrid,sdtime[i][0])
        a5.plot(pitchgrid,sdpitch[i][0])

    a1.set_xlabel('Energy (eV)'), a1.set_ylabel('f(E) (a.u.)'), \
        a1.set_title('Thermal final energy')
    a2.set_xlabel('pitch (1)'), a2.set_ylabel('f(xi) (a.u.)'), \
        a2.set_title('Thermal final pitch')
    a3.set_xlabel('Energy (eV)'), a3.set_ylabel('f(E) (a.u.)'), \
        a3.set_title('Slowing-down energy dist')
    a5.set_xlabel('pitch (1)'), a5.set_ylabel('f(xi) (a.u.)'), \
        a5.set_title('Slowing-down final pitch')

    a4.set_xlabel('Time (s)'), a4.set_ylabel('f(t) (a.u.)'), \
        a4.set_title('Slowing-down time')

    a7.pcolor(vpagrid,vpegrid,sdvdist[0]), a7.set_xlabel('vpa (m/s)'), \
        a7.set_ylabel('vpe (m/s)'), a7.set_title('GO: Slowing down dist')
    a8.pcolor(vpagrid,vpegrid,sdvdist[1]), a8.set_xlabel('vpa (m/s)'), \
        a8.set_ylabel('vpe (m/s)'), a7.set_title('GCF: Slowing down dist')
    a9.pcolor(vpagrid,vpegrid,sdvdist[2]), a9.set_xlabel('vpa (m/s)'), \
        a9.set_ylabel('vpe (m/s)'), a7.set_title('GCA: Slowing down dist')

    l1, = a6.plot(0, 0, label='Analytic')
    l2, = a6.plot(0, 0, label='GO')
    l3, = a6.plot(0, 0, label='GC Fixed')
    l4, = a6.plot(0, 0, label='GC Adaptive')
    handles= [l1, l2, l3, l4]
    labels = ['Analytic', 'GO', 'GC Fixed', 'GC Adaptive']
    a6.legend(handles, labels)
    a6.get_xaxis().set_ticks([])
    a6.get_yaxis().set_ticks([])

    plt.show()

def heaviside(x,x0):
    y = np.zeros(x.shape)
    y[x>=x0] = 1.0
    return y

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
