#!/usr/bin/env python3

"""
Test the evaluation and implementation of the boozer coordinates.

This file initializes, runs, and plots test case for checking that the tools in
a5py correctly construct boozer coordinates, and that these coordinates are
correctly implemented in ASCOT5.

Boozer coordinates have the following properties:

1. Field lines appear as straight lines whose slope equals to the safety factor.

2. The Jacobian has the following dependency h(psi) / B^2.

These properties are tested by implementing an analytical equilibrium,
consturcting the boozer coordinates, and evaluating the field line trajectories
with ASCOT5. Furthermore, it is shown that the helical perturbations appear
correctly.

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
import a5py.ascot5io.mrk_fl    as mrk
import a5py.ascot5io.boozer    as boozer
import a5py.ascot5io.mhd       as mhdmod

import a5py.testascot.helpers as helpers

from a5py.ascotpy.ascotpy import Ascotpy
from a5py.preprocessing.analyticequilibrium import psi0 as psifun
from a5py.preprocessing.generateboozer import generate as genbooz

psi_mult  = 200
R0        = 6.2
z0        = 0
Bphi0     = 5.3

# ITER-like but circular equilibrium
psi_coeff = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                       6.200e-03, -1.205e-03, -3.701e-05,  0,
                       0,          0,          0,          0,         -0.155])

# MHd modes
nmode = 2
mmode = 3

def init():
    """
    Initialize tests

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*             Generate options for field line tracing                     #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 4
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = 1e3
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-1
    odict["ORBITWRITE_NPOINT"]         = 10002

    options.write_hdf5(helpers.testfn, odict, desc="BOOZER")


    #**************************************************************************#
    #*                            Field line markers                           #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 2
    ids    = np.array([1, 2])
    R      = np.array([7.0, 8.0])
    weight = 1       * np.array([1, 1])
    pitch  = 1       * np.array([1, 1])
    time   = 0       * np.array([1, 1])
    phi    = 0       * np.array([1, 1])
    z      = 0       * np.array([1, 1])
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, R, phi, z, pitch, weight, time,
                   desc="BOOZER")

    #**************************************************************************#
    #*                     Construct ITER-like magnetic field                  #
    #*                                                                         #
    #**************************************************************************#
    B_GS.write_hdf5(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                    desc="BOOZER")

    #**************************************************************************#
    #*                     Construct the boozer input                          #
    #*                                                                         #
    #**************************************************************************#
    h5 = ascot5.Ascot(helpers.testfn)
    a5 = Ascotpy(helpers.testfn)
    a5.init(bfield=h5.bfield["BOOZER"].get_qid())
    booz, rhoprof, qprof = genbooz(a5,
                                   rgrid=np.linspace(4, 8.5, 200),
                                   zgrid=np.linspace(-2.5, 2.5, 200),
                                   npsi=200,
                                   nthgeo=360,
                                   nthbzr=360,
                                   raxis=h5.bfield["BOOZER"].read()["raxis"],
                                   zaxis=h5.bfield["BOOZER"].read()["zaxis"],
                                   psi0=h5.bfield["BOOZER"].read()["psi0"],
                                   psi1=h5.bfield["BOOZER"].read()["psi1"])
    a5.free(bfield=True)
    boozer.write_hdf5(helpers.testfn, desc="BOOZER", **booz)

    # Get qprofile value easily by uncommenting this line
    #print(np.interp(np.array([0.21, 0.85]), rhoprof.ravel(), qprof.ravel()))

    #**************************************************************************#
    #* Helical perturbation which is not used in the simulation but is used    #
    #* during the postprocessing                                               #
    #*                                                                         #
    #**************************************************************************#
    mhd = {}
    mhd["nmode"]     = 1
    mhd["nmodes"]    = np.array([nmode])
    mhd["mmodes"]    = np.array([mmode])
    mhd["amplitude"] = np.array([1])
    mhd["omega"]     = np.array([1])

    mhd["npsi"]   = 100
    mhd["psimin"] = 0.0
    mhd["psimax"] = 1

    psigrid = np.linspace(mhd["psimin"], mhd["psimax"], mhd["npsi"])
    alpha = np.exp( -(psigrid-0.85)**2/0.1 )
    phi = alpha*0
    mhd["phi"]   = np.tile(phi, (mhd["nmode"],1)).T
    mhd["alpha"] = np.tile(alpha, (mhd["nmode"],1)).T
    mhdmod.write_hdf5(helpers.testfn, desc="BOOZER", **mhd)

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="BOOZER")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc="BOOZER")
    N0_3D.write_hdf5_dummy(helpers.testfn, desc="BOOZER")

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
                    edens, etemp, idens, itemp, desc="BOOZER")


def run():
    """
    Run tests.
    """
    for test in ["BOOZER"]:
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
    orb = h5["BOOZER"].orbit

    a5 = Ascotpy(helpers.testfn)
    a5.init(bfield=h5.bfield["BOOZER"].get_qid(),
            boozer=h5.boozer["BOOZER"].get_qid(),
            mhd=h5.mhd["BOOZER"].get_qid())
    theta = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                        t=orb["time"], quantity="theta")
    zeta = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                       t=orb["time"], quantity="zeta")
    alpha = a5.evaluate(orb["r"], phi=orb["phi"], z=orb["z"],
                        t=orb["time"]*0, quantity="alpha")
    a5.free(bfield=True, boozer=True, mhd=True)

    ids = orb["id"]
    th = np.linspace(0,2*np.pi,200)

    fig = plt.figure()
    s1 = fig.add_subplot(2,2,1)
    s1.plot(theta[ids==1], np.mod(zeta[ids==1]+2*np.pi, 2*np.pi),
            linestyle="none", marker=".")
    line = 1.67*th
    s1.plot(th[line < 2*np.pi], line[line < 2*np.pi], color="black")

    s2 = fig.add_subplot(2,2,2)
    s2.plot(theta[ids==2], np.mod(zeta[ids==2]+2*np.pi, 2*np.pi),
            linestyle="none", marker=".")
    line = 2.19*th
    s2.plot(th[line < 2*np.pi], line[line < 2*np.pi], color="black")

    s3 = fig.add_subplot(2,2,3)
    s3.plot(np.mod(nmode*zeta[ids==1] - mmode*theta[ids==1], 2*np.pi),
            alpha[ids==1], linestyle="none", marker=".")

    s4 = fig.add_subplot(2,2,4)
    s4.plot(np.mod(nmode*zeta[ids==2] - mmode*theta[ids==2], 2*np.pi),
            alpha[ids==2], linestyle="none", marker=".")

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
