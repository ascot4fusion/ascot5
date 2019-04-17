"""
Test conservation of quantity K = H- omega*P_phi/n for simulations involving
mhd modes. The conservation is tested for each of the orbit integrators:

1. gyro-orbit integrator (VPA) 
2. guiding center fixed step scheme (RK4) 
3. guiding center adaptive step scheme (Cash-Karp)

Tests are done in ITER-like but analytic (B_GS) magnetic field without
collisions. Test particle is a energetic electron (and positron) so tests also
verify that ASCOT5 is valid in relativistic regime.

To init, run and check this test, call this script without any arguments. To
do only one of the above, call this script with an argument "init", "run", or
"check".

File: test_mhd.py
"""
import sys

import numpy                   as np
import scipy.constants         as constants
import matplotlib.pyplot       as plt

import a5py.ascot5io.ascot5    as ascot5
import a5py.ascot5io.orbits    as orbits
import a5py.ascot5io.options   as options
import a5py.ascot5io.B_GS      as B_GS
import a5py.ascot5io.B_2DS     as B_2DS
import a5py.ascot5io.E_TC      as E_TC
import a5py.ascot5io.plasma_1D as P_1D
import a5py.ascot5io.wall_2D   as W_2D
import a5py.ascot5io.N0_3D     as N0_3D
import a5py.ascot5io.mrk_gc    as mrk
import a5py.ascot5io.boozer    as boozer
import a5py.ascot5io.mhd       as mhd
import a5py.preprocessing.eqdsk2input as eqdsk

from a5py.preprocessing.boozermaps import Boozermaps

from a5py.ascotpy import Ascotpy

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
    - MHD-GO tests VPA algorithm used in integrating gyro-orbit motion
    - MHD-GCF tests RK4 used in integrating guiding center motion with fixed
      time-step
    - MHD-GCA tests Cash-Karp used in integrating guiding center motion with
      adaptive time-step

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """

    #**************************************************************************#
    #*                     Generate options for MHD-GO                      #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 1
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-11
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10
    odict["ORBITWRITE_NPOINT"]         = 50002
    odict["ENABLE_MHD"]                = 1

    options.write_hdf5(helpers.testfn, odict, desc="MHD_GO")

    #**************************************************************************#
    #*                     Generate options for MHD-GCF                     #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-10
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8
    odict["ORBITWRITE_NPOINT"]         = 502
    odict["ENABLE_MHD"]                = 1

    options.write_hdf5(helpers.testfn, odict, desc="MHD_GCF")

    #**************************************************************************#
    #*                     Generate options for MHD-GCA                     #
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
    odict["ENDCOND_MAX_SIMTIME"]       = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-8
    odict["ORBITWRITE_NPOINT"]         = 502
    odict["ENABLE_MHD"]                 = 1

    options.write_hdf5(helpers.testfn, odict, desc="MHD_GCA")

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
    anum   = 1       * np.array([1, 0])
    znum   = 1       * np.array([1, 0])
    time   = 0       * np.array([1, 1])
    R      = 0.1     * np.array([1, 1])
    phi    = 90      * np.array([1, 1])
    z      = 0       * np.array([1, 1])
    zeta   = 2       * np.array([1, 1])
    energy = 10e6    * np.array([1, 1])
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="MHD_GO")
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="MHD_GCF")
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="MHD_GCA")

    #**************************************************************************#
    #*                     Construct a field from geqdsk;                      #
    #*                         add MHD modes                                   #
    #**************************************************************************#
    #B_GS.write_hdf5(helpers.testfn, R0,z0, Bphi0, psi_mult, psi_coeff,
#desc="MHD_GO")

    #B_GS.write_hdf5(helpers.testfn, R0,z0, Bphi0, psi_mult, psi_coeff,
#desc="MHD_GCF")


    #B_GS.write_hdf5(helpers.testfn, R0,z0, Bphi0, psi_mult, psi_coeff,
#desc="MHD_GCA")
    B2D = eqdsk.generate("input.geq",cocosin=3,desc="MHD_GO", fn_out=helpers.testfn)

    print(B2D)

    eqdsk.generate("input.geq",cocosin=3,desc="MHD_GCF", fn_out=helpers.testfn)

    eqdsk.generate("input.geq",cocosin=3,desc="MHD_GCA", fn_out=helpers.testfn)
 
    fmhdname = "boozer_maps.out"
    b = Boozermaps(fmhdname)
    b.write_hdf5(helpers.testfn, desc="MHD_GO")
    b.write_hdf5(helpers.testfn, desc="MHD_GCF")
    b.write_hdf5(helpers.testfn, desc="MHD_GCA")

    #create dummy data
    data = {}
    data["nmode"] = 1 
    data["npsi"] = 500
    data["mmodes"] = np.array([10])
    data["nmodes"] = np.array([11])
    data["amplitude"] = np.array([.001])
    data["omega"] = np.array([3.24e+05])
    data["psimin"] = 0.0
    data["psimax"] = 1.0
    psi = np.linspace(data["psimin"],data["psimax"],data["npsi"])
    data["alpha"] = np.zeros((data["npsi"],data["nmode"]))
    for i in range(data["npsi"]):
        data["alpha"][i,0] = np.exp(-50*np.square(psi[i]-0.5))
    data["phi"] = np.zeros((data["npsi"],data["nmode"]))
    for i in range(data["npsi"]):
        data["phi"][i,0] = np.exp(-50*np.square(psi[i]-0.5)) 


    #write that data to file

    mhd.write_hdf5(helpers.testfn, data["nmode"], data["nmodes"],
                   data["mmodes"], data["amplitude"], data["omega"], 
                   data["alpha"], data["phi"],data["npsi"], data["psimin"],
                   data["psimax"], desc="MHD_GO")


    mhd.write_hdf5(helpers.testfn, data["nmode"], data["nmodes"],
                   data["mmodes"], data["amplitude"], data["omega"], 
                   data["alpha"], data["phi"],data["npsi"], data["psimin"],
                   data["psimax"], desc="MHD_GCF")


    mhd.write_hdf5(helpers.testfn, data["nmode"], data["nmodes"],
                   data["mmodes"], data["amplitude"], data["omega"], 
                   data["alpha"], data["phi"],data["npsi"], data["psimin"],
                   data["psimax"], desc="MHD_GCA")
    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="MHD_GO")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="MHD_GCF")
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="MHD_GCA")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    for tname in ["MHD_GO", "MHD_GCF", "MHD_GCA"]:
        W_2D.write_hdf5(helpers.testfn, nwall, Rwall, zwall, desc=tname)
        N0_3D.write_hdf5_dummy(helpers.testfn, desc=tname)

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
                    edens, etemp, idens, itemp, desc="MHD_GO")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="MHD_GCF")
    P_1D.write_hdf5(helpers.testfn, Nrho, Nion, znum, anum, mass, charge, rho,
                    edens, etemp, idens, itemp, desc="MHD_GCA")


def run():
    """
    Run tests.
    """
    for test in ["MHD_GO", "MHD_GCF", "MHD_GCA"]:
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
    apy = Ascotpy(helpers.testfn)
    print("initializing bfield")
    apy.init(bfield=True)
    print("initializing boozer")
    apy.init(boozer=True)
    print("initializing mhd")
    apy.init(mhd=True) 
    
    #must have only one nmode and omega
    nmode = float(a5["MHD_GO"].mhd.read()["nmodes"][0])
    omega = float(a5["MHD_GO"].mhd.read()["omega"][0])
   
    raxis = R0#a5["MHD_GO"].bfield.read()["raxis"]

    f = plt.figure(figsize=(11.9/2.54, 8/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    h1 = f.add_subplot(1,1,1)

    #colors = ["#373e02", "#0a481e", "#03719c", "#0165fc", "#7e1e9c", "#cea2fd"]
    colors = ["b", "dodgerblue", "darkgreen", "forestgreen", "r", "tomato"]

    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for MHD_GO             #
    #*                                                                         #
    #**************************************************************************#
    MHD = {}
    MHD["GO"] = {}
    orb = a5["MHD_GO"]["orbit"].read()

    B = np.sqrt( orb["br"] * orb["br"] + orb["bphi"] * orb["bphi"] +
                 orb["bz"] * orb["bz"] )

    psi = psifun(orb["r"]/raxis, orb["z"]/raxis, psi_coeff[0], psi_coeff[1],
                 psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                 psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                 psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult

    vnorm = np.sqrt( orb["vr"]   * orb["vr"] +
                     orb["vphi"] * orb["vphi"] +
                     orb["vz"]   * orb["vz"] )

    vpar  = ( orb["vr"] * orb["br"] + orb["vphi"] * orb["bphi"] +
              orb["vz"] * orb["bz"] ) / B

    gamma = np.sqrt(1 / ( 1 - vnorm * vnorm / (c * c) ) )

    MHD["GO"]["time"] = orb["time"]
    MHD["GO"]["id"]   = orb["id"]
    MHD["GO"]["r"]    = orb["r"]
    MHD["GO"]["z"]    = orb["z"]
    MHD["GO"]["ekin"] = (gamma - 1) * m_e * c * c
    MHD["GO"]["mu"]   = ( ( m_e * gamma * gamma ) / ( 2 * B ) ) * \
                           ( vnorm * vnorm - vpar * vpar )
    MHD["GO"]["ctor"] = gamma * m_e * orb["r"] * orb["vphi"] + \
                           orb["charge"] * e * psi
    MHD["GO"]["phi"]  = orb["phi"]

   
    id1 = MHD["GO"]["id"] == 1
    id2 = MHD["GO"]["id"] == 2
    print((omega/nmode)*MHD["GO"]["ctor"][id1])
    print(type(omega))
    print(type(nmode))
    print(type(omega/nmode))
    print(apy.evaluate(orb["r"][id1],orb["phi"][id1],orb["z"][id1],orb["time"][id1], "mhd_phi"))
 
    plot_relerr(h1, MHD["GO"]["time"][id1], (MHD["GO"]["ekin"][id1] +
               e*apy.evaluate(orb["r"][id1],orb["phi"][id1],orb["z"][id1], orb["time"][id1], "mhd_phi") 
              -(omega/nmode)*MHD["GO"]["ctor"][id1]), colors[0],label="GO 1")   

    plot_relerr(h1, MHD["GO"]["time"][id2], MHD["GO"]["ekin"][id2] +
               e*apy.evaluate(orb["r"][id2],orb["phi"][id2],orb["z"][id2],
               orb["time"][id2], "mhd_phi") 
               - (omega/nmode)*MHD["GO"]["ctor"][id2], colors[1], label = "GO 2")   

    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for MHD_GCF            #
    #*                                                                         #
    #**************************************************************************#
    MHD["GCF"] = {}
    orb = a5["MHD_GCF"]["orbit"].read()

    B = np.sqrt(np.power(orb["br"],2) + np.power(orb["bphi"],2) +
                np.power(orb["bz"],2))

    psi = psifun(orb["r"]/raxis, orb["z"]/raxis, psi_coeff[0], psi_coeff[1],
                 psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                 psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                 psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult

    # Note that mu is in eV / T
    gamma = np.sqrt( ( 1 + 2 * orb["mu"] * e * B / ( m_e * c * c ) ) /
                     ( 1 - orb["vpar"] * orb["vpar"] / ( c * c ) ) )

    MHD["GCF"]["time"] = orb["time"]
    MHD["GCF"]["id"]   = orb["id"]
    MHD["GCF"]["r"]    = orb["r"]
    MHD["GCF"]["z"]    = orb["z"]
    MHD["GCF"]["ekin"] = (gamma - 1) * m_e * c * c
    MHD["GCF"]["mu"]   = orb["mu"] * e
    MHD["GCF"]["ctor"] = gamma * m_e * orb["r"] * orb["vpar"] + \
                            orb["charge"] * e * psi

    MHD["GCF"]["phi"]  = orb["phi"]
    
    
    id1 = MHD["GCF"]["id"] == 1
    id2 = MHD["GCF"]["id"] == 2

    #plot_relerr(h1, MHD["GCF"]["time"][id1], MHD["GCF"]["ekin"][id1] +
    #           e*apy.evaluate(orb["r"][id1],orb["phi"][id1],orb["z"][id1], orb["time"][id1], "mhd_phi") 
    #           - (omega/nmode)*MHD["GCF"]["ctor"][id1], colors[2],label="GCF 1")  
    
    #plot_relerr(h1, MHD["GCF"]["time"][id2], MHD["GCF"]["ekin"][id2] +
    #           e*apy.evaluate(orb["r"][id2],orb["phi"][id2],orb["z"][id2],
    #           orb["time"][id2], "mhd_phi") 
    #           - (omega/nmode)*MHD["GCF"]["ctor"][id2], colors[3], label="GCF 2")  

    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for MHD_GCA            #
    #*                                                                         #
    #**************************************************************************#
    MHD["GCA"] = {}
    orb = a5["MHD_GCA"]["orbit"].read()

    B = np.sqrt(np.power(orb["br"],2) + np.power(orb["bphi"],2) +
                np.power(orb["bz"],2))

    psi = psifun(orb["r"]/raxis, orb["z"]/raxis, psi_coeff[0], psi_coeff[1],
                 psi_coeff[2], psi_coeff[3], psi_coeff[4], psi_coeff[5],
                 psi_coeff[6], psi_coeff[7], psi_coeff[8], psi_coeff[9],
                 psi_coeff[10], psi_coeff[11], psi_coeff[12]) * psi_mult

    # Note that mu is in eV / T
    gamma = np.sqrt( ( 1 + 2 * orb["mu"] * e * B / ( m_e * c * c ) ) /
                     ( 1 - orb["vpar"] * orb["vpar"] / ( c * c ) ) )

    MHD["GCA"]["time"] = orb["time"]
    MHD["GCA"]["id"]   = orb["id"]
    MHD["GCA"]["r"]    = orb["r"]
    MHD["GCA"]["z"]    = orb["z"]
    MHD["GCA"]["ekin"] = (gamma - 1) * m_e * c * c
    MHD["GCA"]["mu"]   = orb["mu"] * e
    MHD["GCA"]["ctor"] = gamma * m_e * orb["r"] * orb["vpar"] + \
                            orb["charge"] * e * psi

    MHD["GCA"]["phi"] = orb["phi"]


    id1 = MHD["GCA"]["id"] == 1
    id2 = MHD["GCA"]["id"] == 2

    #plot_relerr(h1, MHD["GCA"]["time"][id1], MHD["GCA"]["ekin"][id1] -
    #           e*apy.evaluate(orb["r"][id1],orb["phi"][id1],orb["z"][id1], orb["time"][id1], "mhd_phi") 
    #           -(omega/nmode)*MHD["GCA"]["ctor"][id1], colors[4],label = "GCA 1")  


    #plot_relerr(h1, MHD["GCA"]["time"][id2], MHD["GCA"]["ekin"][id2] +
    #           e*apy.evaluate(orb["r"][id2],orb["phi"][id2],orb["z"][id2],
    #           orb["time"][id2], "mhd_phi") 
    #           - (omega/nmode)*MHD["GCA"]["ctor"][id2], colors[5], label = "GCA 2")  


    #**************************************************************************#
    #*                 Finalize and print and show the figure                  #
    #*                                                                         #
    #**************************************************************************#

    h1.set_xlim(0, 5e-6)
    h1.xaxis.set(ticks=[0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6], ticklabels=[])
    h1.tick_params(axis='y', direction='out')
    h1.tick_params(axis='x', direction='out')
    h1.spines['right'].set_visible(False)
    h1.spines['top'].set_visible(False)
    h1.yaxis.set_ticks_position('left')
    h1.xaxis.set_ticks_position('bottom')
    h1.set(ylabel=r"$\Delta E/E_0$")
    h1.legend()

    plt.savefig("test_mhd.png", dpi=300)
    plt.show()

def plot_relerr(axis, x, y, color,label=None):
    axis.plot(x, y/y[0] - 1, color,label=label)


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
