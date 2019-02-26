"""
Test guiding center transformation.

This file initializes, runs, and plots test case for checking that ASCOT5
performs firstorder guiding center correctly. Since there is no analytical
way to verify this, except for comparing to the very same formulas that ASCOT5
implements, instead we:

1. Initialize a guiding center marker and trace it in GC mode ITER-like
   axisymmetric analytic (B_GS) field without collisions.
2. Trace the very same guiding center with FO mode, so that we get to perform
   guiding center to particle transformation.
3. Repeat step 2. but this time perform particle to guiding center at each time-
   step, and record the guiding center coordinates while continuing the FO
   simulation. This way we can compare the "actual" guiding center trajectory
   to the one computed with GC mode.

Ideally, the actual guiding center trajectory should be close to the GC
trajectory and the former should have almost eliminated the gyroangle
dependency. Do note that the guiding center transformation is not an exact but
done to the first order. Also note that because magnetic moment is not an exact
invariant, the actual guiding center and GC trajectories never match completely.
Test particle is a energetic electron so this test verifies that the
transformation is valid also in the relativstic regime.

This test also shows the difference in first and zeroth order velocity space*
transformation. This is done by taking the calculated orbit from step 2. and
using that to launch additional GC simulations. These simulations are launched
within the first gyroperiod, so ideally the resulting GC trajectories should
be identical.

To init, run and check this test, call this script without any arguments. To
do only one of the above, call this script with an argument "init", "run", or
"check".

*Zeroth order spatial transformation is trivial as then guiding center and
particle coordinates would be identical. Therefore that is not examined here.

File: test_gctransform.py
"""
import sys

import numpy                   as np
import scipy.constants         as constants
import matplotlib.pyplot       as plt
import matplotlib.lines        as mlines

import a5py.ascot5io.ascot5    as ascot5
import a5py.ascot5io.orbits    as orbits
import a5py.ascot5io.options   as options
import a5py.ascot5io.B_GS      as B_GS
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.mrk_gc    as mrk
import a5py.ascot5io.mrk_prt   as prt

import a5py.testascot.helpers as helpers

from a5py.preprocessing.analyticequilibrium import psi0 as psifun

e       = constants.elementary_charge
m_a_AMU = constants.physical_constants["alpha particle mass in u"][0]
m_a     = constants.physical_constants["alpha particle mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

psi_mult  = 200
R0        = 6.2
z0        = 0
Bphi0     = 5.3

# ITER-like equilibrium
psi_coeff = np.array([ 8.629e-02,  3.279e-01,  5.268e-01, -2.366e-01,
                       3.825e-01, -3.573e-01, -1.484e-02,  1.506e-01,
                       7.428e-01, -4.447e-01, -1.084e-01,  1.281e-02, -0.155])

nrep = 10

def init():
    """
    Initialize tests

    This function initializes four test cases:
    - GCTRANSFORM_GC initializes a guiding center marker and simulates it in GC
      mode.
    - GCTRANSFORM_GO initializes a guiding center marker and simulates it in FO
      mode.
    - GCTRANSFORM_GO2GC is same as GO but guiding center transformation is done
      at each
      time-step to get the actual guiding center trajectory.
    - GCTRANSFORM_FIRST initializes particle markers from the orbit calculated
      in GCTRANSFORM_GO. Markers are initialized at different time but all are
      launched within the first gyroperiod. Markers are simulated in GC mode
      with the initial guiding center transformation done to first order.
    - GCTRANSFORM_ZEROTH is same as GCTRANSFORM_FIRST but the guiding center
      transformation is only done to zeroth order in velocity space.

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                     Generate options for GC runs                        #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 1.5e-5
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 2e-10
    odict["ORBITWRITE_MAXPOINTS"]      = 75002

    options.write_hdf5(helpers.testfn, odict, desc="GCTRANSFORM_GC")
    options.write_hdf5(helpers.testfn, odict, desc="GCTRANSFORM_FIRST")

    odict["DISABLE_FIRSTORDER_GCTRANS"] = 1
    options.write_hdf5(helpers.testfn, odict, desc="GCTRANSFORM_ZEROTH")

    #**************************************************************************#
    #*                     Generate options for GO                             #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 1.5e-5
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 2e-10
    odict["ORBITWRITE_MAXPOINTS"]      = 75002

    options.write_hdf5(helpers.testfn, odict, desc="GCTRANSFORM_GO")

    #**************************************************************************#
    #*                     Generate options for GO2GC                          #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["RECORD_GO_AS_GC"]           = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 1.5e-5
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 2e-10
    odict["ORBITWRITE_MAXPOINTS"]      = 75002

    options.write_hdf5(helpers.testfn, odict, desc="GCTRANSFORM_GO2GC")

    #**************************************************************************#
    #*                 Marker input consisting of an alpha particle            #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 1
    ids    = 1
    weight = 1
    pitch  = 0.4
    mass   = m_a_AMU
    charge = 2
    time   = 0
    R      = 7.6
    phi    = 90
    z      = 0
    theta  = 2
    energy = 10e6
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GCTRANSFORM_GC")
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GCTRANSFORM_GO")
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GCTRANSFORM_GO2GC")

    #**************************************************************************#
    #*                     Construct ITER-like magnetic field                  #
    #*                                                                         #
    #**************************************************************************#
    Rmin = 4; Rmax = 8.5; nR = 120; zmin = -4; zmax = 4; nz = 200;
    B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_GC")
    B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_GO")
    B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_GO2GC")
    B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_ZEROTH")
    B_GS.write_hdf5_B_2D(helpers.testfn, R0, z0, Bphi0, psi_mult, psi_coeff,
                         Rmin, Rmax, nR, zmin, zmax, nz,
                         desc="GCTRANSFORM_FIRST")

    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="GCTRANSFORM_GC")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="GCTRANSFORM_GO")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="GCTRANSFORM_GO2GC")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="GCTRANSFORM_ZEROTH")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="GCTRANSFORM_FIRST")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_GC")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_GO")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_GO2GC")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_ZEROTH")
    W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall,
                    desc="GCTRANSFORM_FIRST")

    helpers.write_N0_3D_dummy(helpers.testfn, desc="GCTRANSFORM_GC")
    helpers.write_N0_3D_dummy(helpers.testfn, desc="GCTRANSFORM_GO")
    helpers.write_N0_3D_dummy(helpers.testfn, desc="GCTRANSFORM_GO2GC")
    helpers.write_N0_3D_dummy(helpers.testfn, desc="GCTRANSFORM_ZEROTH")
    helpers.write_N0_3D_dummy(helpers.testfn, desc="GCTRANSFORM_FIRST")

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
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_GC")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_GO")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_GO2GC")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_ZEROTH")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GCTRANSFORM_FIRST")

def run():
    """
    Run tests.

    This first runs the GCTRANSFORM_GC, GCTRANSFORM_GO, and GCTRANSFORM_GO2GC
    cases, then uses results of GCTRANSFORM_GO to generate markers for
    GCTRANSFORM_ZEROTH and GCTRANSFORM_FIRST which are then run.
    """
    for test in ["GCTRANSFORM_GC", "GCTRANSFORM_GO", "GCTRANSFORM_GO2GC"]:
        helpers.set_and_run(test)

    a5 = ascot5.Ascot(helpers.testfn)

    dt = 20
    Nmrk   = nrep
    ids    = np.linspace(1, Nmrk, Nmrk)
    weight = 1       * np.ones(ids.shape)
    mass   = m_a_AMU * np.ones(ids.shape)
    charge = 2       * np.ones(ids.shape)
    time   = a5["GCTRANSFORM_GO"]["orbits"]["time"][0:Nmrk*dt:dt]
    R      = a5["GCTRANSFORM_GO"]["orbits"]["R"][0:Nmrk*dt:dt]
    phi    = a5["GCTRANSFORM_GO"]["orbits"]["phi"][0:Nmrk*dt:dt]
    z      = a5["GCTRANSFORM_GO"]["orbits"]["z"][0:Nmrk*dt:dt]
    vR     = a5["GCTRANSFORM_GO"]["orbits"]["v_R"][0:Nmrk*dt:dt]
    vphi   = a5["GCTRANSFORM_GO"]["orbits"]["v_phi"][0:Nmrk*dt:dt]
    vz     = a5["GCTRANSFORM_GO"]["orbits"]["v_z"][0:Nmrk*dt:dt]
    prt.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, vR, vphi, vz,
                   weight, time, desc="GCTRANSFORM_ZEROTH")
    prt.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, vR, vphi, vz,
                   weight, time, desc="GCTRANSFORM_FIRST")

    for test in ["GCTRANSFORM_ZEROTH", "GCTRANSFORM_FIRST"]:
        helpers.set_and_run(test)

def check():
    """
    Plot the results of these tests.

    This function makes four plots.
    - One that shows difference in magnetic moment with respect to time between
      cases GCTRANSFORM_GC and GCTRANSFORM_GO, and GCTRANSFORM_GC and
      GCTRANSFORM_GO2GC.
    - The second one is the same but for the parallel velocity.
    - Third plot shows the orbits of these three cases on a Rz plane with focus
      on the banana turning point.
    - Fourth one is the same as previous plot but showing cases GCTRANSFORM_GC,
      GCTRANSFORM_ZEROTH and GCTRANSFORM_FIRST. The latter two have orbits from
      multiple markers but these are plotted with same color.
    """
    a5 = ascot5.Ascot(helpers.testfn)

    f = plt.figure(figsize=(11.9/2.54, 8/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    h1 = f.add_subplot(1,4,1)
    h1.set_position([0.12, 0.56, 0.26, 0.4], which='both')

    h2 = f.add_subplot(1,4,2)
    h2.set_position([0.12, 0.15, 0.26, 0.4], which='both')

    h3 = f.add_subplot(1,4,3)
    h3.set_position([0.5, 0.15, 0.24, 0.8], which='both')

    h4 = f.add_subplot(1,4,4)
    h4.set_position([0.74, 0.15, 0.24, 0.8], which='both')

    # Make plots, scale the plotted quantities.
    # For some reason GCTRANSFORM_GO2GC GCTRANSFORM_GC are not equal in length
    # and we need to have [:-1]?
    h1.plot(a5["GCTRANSFORM_GO"]["orbits"]["time"]*1e6,
            ( a5["GCTRANSFORM_GO"]["orbits"]["mu"]
              - a5["GCTRANSFORM_GC"]["orbits"]["mu"] )/e / 1e4 )
    h1.plot(a5["GCTRANSFORM_GO2GC"]["orbits"]["time"]*1e6,
            ( a5["GCTRANSFORM_GO2GC"]["orbits"]["mu"]
              - a5["GCTRANSFORM_GC"]["orbits"]["mu"][:-1] )/e / 1e4 )

    h2.plot(a5["GCTRANSFORM_GO"]["orbits"]["time"]*1e6,
            ( a5["GCTRANSFORM_GO"]["orbits"]["vpar"]
              - a5["GCTRANSFORM_GC"]["orbits"]["vpar"] ) / 1e5 )
    h2.plot(a5["GCTRANSFORM_GO2GC"]["orbits"]["time"]*1e6,
            ( a5["GCTRANSFORM_GO2GC"]["orbits"]["vpar"]
              - a5["GCTRANSFORM_GC"]["orbits"]["vpar"][:-1] ) / 1e5 )

    h3.plot(a5["GCTRANSFORM_GO"]["orbits"]["R"],
            a5["GCTRANSFORM_GO"]["orbits"]["z"])
    h3.plot(a5["GCTRANSFORM_GO2GC"]["orbits"]["R"],
            a5["GCTRANSFORM_GO2GC"]["orbits"]["z"])
    h3.plot(a5["GCTRANSFORM_GC"]["orbits"]["R"],
            a5["GCTRANSFORM_GC"]["orbits"]["z"])

    for i in range(0, nrep):
        id0 = a5["GCTRANSFORM_ZEROTH"]["orbits"]["id"]
        id1 = a5["GCTRANSFORM_FIRST"]["orbits"]["id"]
        h4.plot(a5["GCTRANSFORM_ZEROTH"]["orbits"]["R"][id0==i+1],
                a5["GCTRANSFORM_ZEROTH"]["orbits"]["z"][id0==i+1], 'darkorchid')
        h4.plot(a5["GCTRANSFORM_FIRST"]["orbits"]["R"][id1==i+1],
                a5["GCTRANSFORM_FIRST"]["orbits"]["z"][id1==i+1], 'darkgreen')

    h4.plot(a5["GCTRANSFORM_GO2GC"]["orbits"]["R"],
            a5["GCTRANSFORM_GO2GC"]["orbits"]["z"], 'red')

    #**************************************************************************#
    #*                 Finalize and print and show the figure                  #
    #*                                                                         #
    #**************************************************************************#

    h1.set_xlim(0, 15)
    h1.set_ylim(-6, 6)
    h1.spines['right'].set_visible(False)
    h1.spines['top'].set_visible(False)
    h1.spines['bottom'].set_visible(False)
    h1.yaxis.set_ticks_position('left')
    h1.xaxis.set_ticks_position('bottom')
    h1.tick_params(axis='y', direction='out')
    h1.tick_params(axis='x', direction='inout')
    h1.xaxis.set(ticks=[0, 5, 10, 15], ticklabels=[])
    h1.yaxis.set(ticks=[-6, -4, -2, 0, 2, 4, 6],
                 ticklabels=[-6, '', '', 0, '', '', 6])
    h1.set(ylabel=r"$\Delta \mu$ [$10^{4}$ eV/T]")

    h2.set_xlim(0, 15)
    h2.set_ylim(-2, 6)
    h2.spines['right'].set_visible(False)
    h2.spines['top'].set_visible(False)
    h2.yaxis.set_ticks_position('left')
    h2.xaxis.set_ticks_position('bottom')
    h2.tick_params(axis='y', direction='out')
    h2.tick_params(axis='x', direction='out')
    h2.xaxis.set(ticks=[0, 5, 10, 15])
    h2.yaxis.set(ticks=[-2, 0, 2, 4, 6],
                 ticklabels=[-2, 0, 2, 4, ''])
    h2.set(ylabel=r"$\Delta v_\parallel$ [$10^{5}$ m/s]",
           xlabel=r"Time [$10^{-6}$ s]")

    h3.axis('scaled')
    h3.set_xlim(6.18, 6.68)
    h3.set_ylim(1.6, 2.6)
    h3.tick_params(axis='y', direction='out')
    h3.tick_params(axis='x', direction='out')
    h3.xaxis.set(ticks=[6.2, 6.4, 6.6])
    h3.yaxis.set(ticks=[1.7, 1.9, 2.1, 2.3, 2.5],
                 ticklabels=[1.7, '', 2.1, '', 2.5])
    h3.set(ylabel=r"$z$ [m]")
    h3.set_xlabel(r"$R$ [m]", position=[1,0,1,1]);

    h4.axis('scaled')
    h4.set_xlim(6.18, 6.68)
    h4.set_ylim(1.6, 2.6)
    h4.tick_params(axis='y', direction='inout')
    h4.tick_params(axis='x', direction='out')
    h4.xaxis.set(ticks=[6.2, 6.4, 6.6])
    h4.yaxis.set(ticks=[1.7, 1.9, 2.1, 2.3, 2.5], ticklabels=[])

    l1 = mlines.Line2D([], [], color='r', linestyle='-', label='GC',    axes=h3)
    l2 = mlines.Line2D([], [], color='b', linestyle='-', label='GO',    axes=h3)
    l3 = mlines.Line2D([], [], color='g', linestyle='-', label='GO2GC', axes=h3)
    h3.legend(handles=[l1, l2, l3], loc='lower left', frameon=False,
              fontsize=9)

    l1 = mlines.Line2D([], [], color='r', linestyle='-', label='GC',    axes=h4)
    l2 = mlines.Line2D([], [], color='darkorchid', linestyle='-',
                       label=r'0$^\circ$ order GC', axes=h4)
    l3 = mlines.Line2D([], [], color='darkgreen', linestyle='-',
                       label=r'1$^\circ$ order GC', axes=h4)
    h4.legend(handles=[l1, l2, l3], loc='lower left', frameon=False,
              fontsize=9)

    plt.savefig("test_gctransform.png", dpi=300)
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
