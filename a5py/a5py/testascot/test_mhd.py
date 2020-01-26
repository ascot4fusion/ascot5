#!/usr/bin/env python3

"""
Test the implementation of the MHD module.

This file initializes, runs, and plots test case for checking that the helical
perturbations introduced by the MHD module conserve the quantity

H = K - omega*P/n,

where K is particle hamiltonian, P canonical angular toroidal momentum, omega
mode frequency and n mode toroidal number.

This is validated for all particle simulation modes.

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
import a5py.ascot5io.boozer    as boozer
import a5py.ascot5io.mhd       as mhdmod

import a5py.testascot.helpers as helpers

from a5py.ascotpy.ascotpy import Ascotpy
from a5py.preprocessing.analyticequilibrium import psi0 as psifun
from a5py.preprocessing.generateboozer import generate as genbooz

e       = constants.elementary_charge
m_p_AMU = constants.physical_constants["proton mass in u"][0]
m_p     = constants.physical_constants["proton mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

psi_mult  = 200
R0        = 6.2
z0        = 0
Bphi0     = 5.3

# ITER-like but circular equilibrium
psi_coeff = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                       6.200e-03, -1.205e-03, -3.701e-05,  0,
                       0,          0,          0,          0,         -0.155])

# MHD mode specifications
nmode = 2
mmode = 3
modefreq = 3e5
amplitude = 1e-2

def init():
    """
    Initialize tests

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*             Generate options for GCF                                    #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-9
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = 1e-3
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_MHD"]                = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-7
    odict["ORBITWRITE_NPOINT"]         = 1e4

    options.write_hdf5(helpers.testfn, odict, desc="MHD_GCF")


    #**************************************************************************#
    #*                            Field line markers                           #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 2
    ids    = np.array([1, 2])
    R      = np.array([7.0, 8.0])
    pitch  = np.array([0.4, 0.9])
    zeta   = 2       * np.array([1, 1])
    energy = 1e5     * np.array([1, 1])
    mass   = m_p_AMU * np.array([1, 1])
    charge = 1       * np.array([1, 1])
    anum   = 1       * np.array([1, 0])
    znum   = 1       * np.array([1, 0])
    weight = 1       * np.array([1, 1])
    time   = 0       * np.array([1, 1])
    phi    = 0       * np.array([1, 1])
    z      = 0       * np.array([1, 1])
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="MHD_GCF")

    #**************************************************************************#
    #*                     Construct ITER-like magnetic field                  #
    #*                                                                         #
    #**************************************************************************#
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="MHD_GCF")

    #**************************************************************************#
    #*                     Construct the boozer input                          #
    #*                                                                         #
    #**************************************************************************#
    h5 = ascot5.Ascot(helpers.testfn)
    a5 = Ascotpy(helpers.testfn)
    a5.init(bfield=h5.bfield["MHD_GCF"].get_qid())
    booz, rhoprof, qprof = genbooz(a5,
                                   rgrid=np.linspace(4, 8.5, 200),
                                   zgrid=np.linspace(-2.5, 2.5, 200),
                                   npsi=200,
                                   nthgeo=360,
                                   nthbzr=360,
                                   raxis=h5.bfield["MHD_GCF"].read()["raxis"],
                                   zaxis=h5.bfield["MHD_GCF"].read()["zaxis"],
                                   psi0=h5.bfield["MHD_GCF"].read()["psi0"],
                                   psi1=h5.bfield["MHD_GCF"].read()["psi1"])
    a5.free(bfield=True)
    boozer.write_hdf5(helpers.testfn, desc="MHD_GCF", **booz)

    #**************************************************************************#
    #* Helical perturbation which is not used in the simulation but is used    #
    #* during the postprocessing                                               #
    #*                                                                         #
    #**************************************************************************#
    mhd = {}
    mhd["nmode"]     = 1
    mhd["nmodes"]    = np.array([nmode])
    mhd["mmodes"]    = np.array([mmode])
    mhd["amplitude"] = np.array([amplitude])
    mhd["omega"]     = np.array([modefreq])

    mhd["npsi"]   = 100
    mhd["psimin"] = 0.0
    mhd["psimax"] = 1

    psigrid = np.linspace(mhd["psimin"], mhd["psimax"], mhd["npsi"])
    alpha = np.exp( -(psigrid-0.85)**2/0.1 )
    phi = alpha*0
    mhd["phi"]   = np.tile(phi, (mhd["nmode"],1)).T
    mhd["alpha"] = np.tile(alpha, (mhd["nmode"],1)).T
    mhdmod.write_hdf5(helpers.testfn, desc="MHD_GCF", **mhd)

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="MHD_GCF")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="MHD_GCF")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="MHD_GCF")

    Nrho   = 3
    Nion   = 1
    znum   = np.array([1])
    anum   = np.array([1])
    mass   = np.array([1])
    charge = np.array([1])
    rho    = np.array([0, 0.5, 100])
    edens  = 1e20 * np.ones(rho.shape)
    etemp  = 1e3  * np.ones(rho.shape)
    idens  = 1e20 * np.ones((rho.size, Nion))
    itemp  = 1e3  * np.ones(rho.shape)
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="MHD_GCF")


def run():
    """
    Run tests.
    """
    for test in ["MHD_GCF"]:
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
    h5 = ascot5.Ascot(helpers.testfn)
    orb = h5["MHD_GCF"].orbit

    B0 = np.sqrt(np.power(orb["br"],2) + np.power(orb["bphi"],2) +
                np.power(orb["bz"],2))

    a5 = Ascotpy(helpers.testfn)
    a5.init(bfield=h5.bfield["MHD_GCF"].get_qid(),
            boozer=h5.boozer["MHD_GCF"].get_qid(),
            mhd=h5.mhd["MHD_GCF"].get_qid())
    psi = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                      t=orb["time"], quantity="psi")
    b0 = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                     t=orb["time"], quantity="bnorm")
    db = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                     t=orb["time"], quantity="db/b")
    dbphi = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                     t=orb["time"], quantity="mhd_bphi")
    a5.free(bfield=True, boozer=True, mhd=True)

    B = b0*(1+db)

    #gamma = np.sqrt( ( 1 + 2 * orb["mu"] * e * B / ( m_p * c * c ) ) /
    #                 ( 1 - orb["vpar"] * orb["vpar"] / ( c * c ) ) )

    #ekin = (gamma - 1) * m_p * c * c

    ekin = orb["mu"]*b0 + 0.5*m_p*orb["vpar"]**2

    #ctor = gamma * m_p * orb["r"] * orb["vpar"] * orb["bphi"] / B + \
    #       orb["charge"] * e * psi
    ctor = m_p * orb["r"] * orb["vpar"] * orb["bphi"] / b0 + \
           orb["charge"] * psi

    ids = orb["id"]
    H = ekin - modefreq*ctor/nmode

    fig = plt.figure()
    plt.plot(H[ids==1], color="black")
    plt.plot(ekin[ids==1], color="red")
    plt.plot(modefreq*ctor[ids==1]/nmode, color="blue")
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
