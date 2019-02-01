"""
Test elementary results such as gyromotion and drifts

This file initializes, runs, and plots test case for checking that ASCOT5
reproduces:

1. correct gyroradius and gyrofrequency in uniform magnetic field
2. correct E X B drift in uniform electromagnetic field
3. correct gradB drift in a magnetic field with constant gradient

Test 1 is done with GO mode and 2 and 3 with both GO and GC modes, latter
using the fixed step scheme. Tests are done in Cartesian electromagnetic
field (B_TC and E_TC) and without collisions. Test particle is a energetic
electron (and positron) so tests also verify that ASCOT5 is valid in
relativistic regime.

To init, run and check this test, call this script without any arguments. To
do only one of the above, call this script with an argument "init", "run", or
"check".

File: test_elementary.py
"""

import sys

import numpy                   as np
import scipy.constants         as constants
import matplotlib.pyplot       as plt

import a5py.ascot5io.ascot5    as ascot5
import a5py.ascot5io.orbits    as orbits
import a5py.ascot5io.options   as options
import a5py.ascot5io.B_TC      as B_TC
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.N0_3D     as N0_3D
import a5py.ascot5io.mrk_gc    as mrk

sys.path.insert(0, '../')
sys.path.insert(0, '.')
import opt
import test_ascot

e       = constants.elementary_charge
m_e_AMU = constants.physical_constants["electron mass in u"][0]
m_e     = constants.physical_constants["electron mass"][0]
c       = constants.physical_constants["speed of light in vacuum"][0]

def init():
    """
    Initialize tests

    This function initializes five test cases:
    - GYROMOTION tests that gyroradius and gyrofrequency are correct
    - EXB-GO tests that ExB-drift is correct for GO
    - EXB-GC tests that ExB-drift is correct for GC
    - GRADB-GO tests that gradB-drift is correct for GO
    - GRADB-GC tests that gradB-drift is correct for GC

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                     Generate options for GYROMOTION                     #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-11
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 2e-9
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-11
    odict["ORBITWRITE_MAXPOINTS"]      = 202

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="GYROMOTION")

    #**************************************************************************#
    #*              Generate options for EXB-GO and GRAD-GO                    #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 1e-7
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-11
    odict["ORBITWRITE_MAXPOINTS"]      = 10002

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="EXB_GO")
    options.write_hdf5(test_ascot.testfn, odict, desc="GRADB_GO")

    #**************************************************************************#
    #*       Generate options for EXB-GC and GRADB-GC                          #
    #*                                                                         #
    #**************************************************************************#
    odict = opt.generateopt()
    test_ascot.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-9
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIM_TIME"]      = 1e-7
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-9
    odict["ORBITWRITE_MAXPOINTS"]      = 102

    opt.settypes(odict)
    options.write_hdf5(test_ascot.testfn, odict, desc="EXB_GC")
    options.write_hdf5(test_ascot.testfn, odict, desc="GRADB_GC")

    #**************************************************************************#
    #*                      Markers are same for all tests                     #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 2
    ids    = np.array([1, 2])
    weight = np.array([1, 1])
    mass   = m_e_AMU * np.array([1, 1])
    charge = 1       * np.array([1,-1])
    time   = 0       * np.array([1, 2])
    R      = 5       * np.array([1, 1])
    phi    = 90      * np.array([1, 1])
    z      = 0       * np.array([1, 1])
    theta  = 0       * np.array([1, 1])
    energy = 100e6   * np.array([1, 1])
    pitch  = 0.5     * np.array([1, 1])
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GYROMOTION")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="EXB_GO")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="EXB_GC")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GRADB_GO")
    mrk.write_hdf5(test_ascot.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, theta,
                   weight, time, desc="GRADB_GC")

    #**************************************************************************#
    #*             Magnetic and electric fields for GYROMOTION                 #
    #*                                                                         #
    #**************************************************************************#
    Bxyz   = np.array([5, 0, 0])
    gradB  = np.array([0,0,0,0,0,0,0,0,0])
    rhoval = 1.5
    B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval, desc="GYROMOTION")

    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="GYROMOTION")

    #**************************************************************************#
    #*            Magnetic and electric fields for EXB-GO and EXB-GC           #
    #*                                                                         #
    #**************************************************************************#
    Bxyz   = np.array([5, 0, 0])
    gradB  = np.array([0,0,0,0,0,0,0,0,0])
    rhoval = 1.5
    B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval, desc="EXB_GO")
    B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval, desc="EXB_GC")

    Exyz   = np.array([0, 1e6, 0])
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="EXB_GO")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="EXB_GC")

    #**************************************************************************#
    #*           Magnetic and electric fields for GRADB-GO and GRADB-GC        #
    #*                                                                         #
    #**************************************************************************#
    Bxyz   = np.array([5, 0, 0])
    gradB  = np.array([0,0,0.1,0,0,0,0,0,0])
    rhoval = 1.5
    B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval, desc="GRADB_GO")
    B_TC.write_hdf5(test_ascot.testfn, Bxyz, gradB, rhoval, desc="GRADB_GC")

    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="GRADB_GO")
    E_TC.write_hdf5(test_ascot.testfn, Exyz, desc="GRADB_GC")

    #**************************************************************************#
    #*              Other inputs are trivial and same for all tests            #
    #*                                                                         #
    #**************************************************************************#
    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="GYROMOTION")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="EXB_GO")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="EXB_GC")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="GRADB_GO")
    W_2D.write_hdf5(test_ascot.testfn, nwall, Rwall, zwall, desc="GRADB_GC")

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
                     N0pmin, N0pmax, N0np, N0dens, desc="GYROMOTION")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="EXB_GO")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="EXB_GC")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="GRADB_GO")
    N0_3D.write_hdf5(test_ascot.testfn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np, N0dens, desc="GRADB_GC")

    Nrho  = 3
    Nion  = 1
    znum  = np.array([1])
    anum  = np.array([1])
    rho   = np.array([0, 0.5, 100])
    edens = 1e20 * np.ones(rho.shape)
    etemp = 1e3  * np.ones(rho.shape)
    idens = 1e20 * np.ones((rho.size, Nion))
    itemp = 1e3  * np.ones(rho.shape)
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GYROMOTION")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="EXB_GO")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="EXB_GC")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GRADB_GO")
    P_1D.write_hdf5(test_ascot.testfn, Nrho, Nion, znum, anum, rho,
                    edens, etemp, idens, itemp, desc="GRADB_GC")

def run():
    """
    Run tests.
    """
    for test in ["GYROMOTION", "EXB_GO", "EXB_GC", "GRADB_GO", "GRADB_GC"]:
        test_ascot.set_and_run(test)

def check():
    """
    Plot the results of these tests.

    This function makes three plots.
    - One that shows gyromotion viewed from above (to chect test GYROMOTION)
    - Two that shows gyromotion viewed from side (to show drifts in EXB-GO,
      EXB_GC, GRADB-GO, and GRADB-GC)
    - Also the analytical an numerical values are shown
    """
    a5 = ascot5.Ascot(test_ascot.testfn)

    GYROMOTION = {}
    EXB        = {}
    GRADB      = {}
    GYROMOTION["GO"] = a5["GYROMOTION"]["orbits"].read()
    EXB["GO"]        = a5["EXB_GO"]["orbits"].read()
    EXB["GC"]        = a5["EXB_GC"]["orbits"].read()
    GRADB["GO"]      = a5["GRADB_GO"]["orbits"].read()
    GRADB["GC"]      = a5["GRADB_GC"]["orbits"].read()

    # Electron energy in Joules
    E = 100e6 * e
    # Lorentz factor
    gamma = 1 + E / ( m_e * c * c )
    # Pitch, magnetic field strength and gradient, electric field, and velocity
    xi    = 0.5
    B     = 5
    gradB = 0.1
    E     = 1e6
    v     = np.sqrt(1.0 - 1.0 / ( gamma * gamma ) ) * c

    f = plt.figure(figsize=(11.9/2.54, 8/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    h1 = f.add_subplot(1,3,1)
    h1.set_position([0.12, 0.27, 0.26, 1], which='both')

    h2 = f.add_subplot(1,3,2)
    h2.set_position([0.56, 0.30, 0.42, 1], which='both')

    h3 = f.add_subplot(1,3,3)
    h3.set_position([0.56, -0.13, 0.42, 1], which='both')

    #**************************************************************************#
    #*                           Check GYROMOTION                              #
    #*                                                                         #
    #**************************************************************************#

    # Analytical values
    rhog   = gamma * np.sqrt(1 - 0.5 * 0.5) * m_e * v / (B * e)
    omegag = e * B / ( gamma * m_e)

    # Numerical values
    ang  = GYROMOTION["GO"]["phi"] * np.pi / 180
    igo  = GYROMOTION["GO"]["id"]
    x    = GYROMOTION["GO"]["R"] * np.sin(ang)
    y    = GYROMOTION["GO"]["z"]
    time = GYROMOTION["GO"]["time"][igo==1]

    rho   = np.max( np.sqrt((x-5)*(x-5) + y*y) )
    zero_crossings = np.where(np.diff(np.sign(x[igo==1])))[0].size
    omega = zero_crossings * np.pi * 2 / time[-1]

    # Plot
    h1.plot(x[igo==1] - 5, y[igo==1])

    #**************************************************************************#
    #*                           Check EXB                                     #
    #*                                                                         #
    #**************************************************************************#

    # Analytical values
    v_ExB = E * B / (B*B)

    # Numerical values
    ang = EXB["GO"]["phi"] * np.pi / 180
    igo = EXB["GO"]["id"]
    xgo = EXB["GO"]["R"] * np.sin(ang)
    ygo = EXB["GO"]["z"]

    igo0  = a5["EXB_GO"]["inistate"]["id"]
    time  = a5["EXB_GO"]["endstate"]["time"]
    ang   = a5["EXB_GO"]["inistate"]["phi"] * np.pi / 180
    xgo0  = a5["EXB_GO"]["inistate"]["R"] * np.sin(ang)
    ygo0  = a5["EXB_GO"]["inistate"]["z"]
    ang   = a5["EXB_GO"]["endstate"]["phi"] * np.pi / 180
    xgo1  = a5["EXB_GO"]["endstate"]["R"] * np.sin(ang)
    ygo1  = a5["EXB_GO"]["endstate"]["z"]

    ang = EXB["GC"]["phi"] * np.pi / 180
    igc = EXB["GC"]["id"]
    xgc = EXB["GC"]["R"] * np.sin(ang)
    ygc = EXB["GC"]["z"]

    vgo1_ExB = ((ygo1[igo0==1] - ygo0[igo0==1]) / (time[igo0==1]))[0]

    y          = ygc[igc==1]
    t          = EXB["GC"]["time"][igc==1]
    vgc1_ExB = (y[-1] - y[0]) / (t[-1] - t[0])

    # Plot
    h2.plot(xgo[igo==1] - 0.07 - 5, ygo[igo==1])
    h2.plot(xgo[igo==2] + 0.07 - 5, ygo[igo==2])
    h2.plot(xgc[igc==1] - 0.07 - 5, ygc[igc==1], color="red")
    h2.plot(xgc[igc==2] + 0.07 - 5, ygc[igc==2], color="red")

    #**************************************************************************#
    #*                           Check GRADB                                   #
    #*                                                                         #
    #**************************************************************************#

    # Analytical values
    v_gradB = ( gamma * m_e * (1.0 - xi*xi) * v * v / ( 2 * e * B) ) \
              * gradB * B / (B*B)

    # Numerical values
    ang = GRADB["GO"]["phi"] * np.pi / 180
    igo = GRADB["GO"]["id"]
    xgo = GRADB["GO"]["R"] * np.sin(ang)
    ygo = GRADB["GO"]["z"]

    igo0  = a5["GRADB_GO"]["inistate"]["id"]
    time  = a5["GRADB_GO"]["endstate"]["time"]
    ang   = a5["GRADB_GO"]["inistate"]["phi"] * np.pi / 180
    xgo0  = a5["GRADB_GO"]["inistate"]["R"] * np.sin(ang)
    zgo0  = a5["GRADB_GO"]["inistate"]["z"]
    ang   = a5["GRADB_GO"]["endstate"]["phi"] * np.pi / 180
    xgo0  = a5["GRADB_GO"]["endstate"]["R"] * np.sin(ang)
    ygo0  = a5["GRADB_GO"]["endstate"]["z"]

    ang = GRADB["GC"]["phi"] * np.pi / 180
    igc = GRADB["GC"]["id"]
    xgc = GRADB["GC"]["R"] * np.sin(ang)
    ygc = GRADB["GC"]["z"]

    vgo1_gradB = ((xgo0[igo0==1] - xgo1[igo0==1]) / (time[igo0==1]))[0]

    x          = xgc[igc==1]
    t          = GRADB["GC"]["time"][igc==1]
    vgc1_gradB = (x[-1] - x[0]) / (t[-1] - t[0])

    # Plot
    h3.plot(ygo[igo==1] - 0.07, xgo[igo==1] - 5)
    h3.plot(ygo[igo==2] + 0.07, xgo[igo==2] - 5)
    h3.plot(ygc[igc==1] - 0.07, xgc[igc==1] - 5, color="red")
    h3.plot(ygc[igc==2] + 0.07, xgc[igc==2] - 5, color="red")

    # Print analytical values
    text1  = r"$\rho_{g}$ = %2.3f cm" % (rhog*100)
    text1 += r" (%2.3f cm)" % (rho*100)
    text1 += "\n"+ r"$\omega_{g}$ = " + latex_float(omegag) + r" Hz "
    text1 += r"(" + latex_float(omegag) + r" Hz)"
    text1 += "\n" + r"$v_{E \times B}$ = " + latex_float(v_ExB) + r" m/s"
    text1 += "\n" + r"(" + latex_float(-vgo1_ExB) + r" m/s, "
    text1 += latex_float(-vgc1_ExB) + r" m/s)" + "\n"
    text1 += r"$v_{\nabla B}$ = " + latex_float(v_gradB) + r" m/s"
    text1 += "\n" + r"(" + latex_float(-vgo1_gradB) + r" m/s, "
    text1 += latex_float(-vgc1_gradB) + r" m/s)" + "\n"
    h1.text(-0.1, -0.24, text1, fontsize=9)

    #**************************************************************************#
    #*                           Finalize                                      #
    #*                                                                         #
    #**************************************************************************#

    # Decorate all plots
    h1.axis('scaled')
    h1.xaxis.set(ticks=[-0.06, 0, 0.06], ticklabels=[-6, 0, 6])
    h1.yaxis.set(ticks=[-0.06, 0, 0.06], ticklabels=[-6, 0, 6])
    h1.tick_params(axis='y', direction='out')
    h1.tick_params(axis='x', direction='out')
    h1.set(xlabel="$x$ [cm]", ylabel="$y$ [cm]")

    h2.axis('scaled')
    h2.set(xlim=[-0.15, 0.15], ylim=[-0.08, 0.08])
    h2.xaxis.set(ticklabels=[])
    h2.yaxis.set(ticks=[-0.08, 0, 0.08], ticklabels=[-8, 0, 8])
    h2.tick_params(axis='y', direction='out')
    h2.tick_params(axis='x', direction='out')
    h2.set(ylabel="$y$ [cm]")

    h3.axis('scaled')
    h3.set(xlim=[-0.15, 0.15], ylim=[-0.08, 0.08])
    h3.xaxis.set(ticks=[-0.15, 0, 0.15], ticklabels=[-15, 0, 15])
    h3.yaxis.set(ticks=[-0.08, 0, 0.08], ticklabels=[-8, 0, 8])
    h3.tick_params(axis='y', direction='out')
    h3.tick_params(axis='x', direction='out')
    h3.set(xlabel="$x$ [cm]", ylabel="$y$ [cm]")

    plt.savefig("test_elementary.png", dpi=300)
    plt.show()

def latex_float(f):
    float_str = "{0:.4g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:
        return float_str

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
