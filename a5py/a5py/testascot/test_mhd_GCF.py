"""
Test conservation of quantity K = H- omega*P_phi/n for simulations involving
mhd modes. The conservation is tested for just one orbit integrator:

1. guiding center fixed step scheme (RK4)

Tests are done using a boozer_maps.out file produced by ASCOT 4 from a .geq 
file. That same .geq file is then read in to create the background equilibrium.

The mode is a dummy mode with frequency, mode numbers, and radial profile
chosen by the tester.

Test particles are an energetic electron and a not-so energetic alpha particle,
 so tests also verify that ASCOT5 is valid in relativistic regime.

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
m_alpha = constants.physical_constants["proton mass"][0]*4 

def init():
    """
    Initialize tests

    This function initializes one test cases:
    - MHD-GCF tests RK4 used in integrating guiding center motion with fixed
      time-step

    Input fields contain the test case name (to which the input corresponds to)
    as a description.
    """


    #**************************************************************************#
    #*                     Generate options for MHD-GCF                     #
    #*                                                                         #
    #**************************************************************************#
    odict = options.generateopt()
    helpers.clean_opt(odict)

    odict["SIM_MODE"]                  = 2
    odict["FIXEDSTEP_USE_USERDEFINED"] = 1
    odict["FIXEDSTEP_USERDEFINED"]     = 1e-12
    odict["ENDCOND_SIMTIMELIM"]        = 1
    odict["ENDCOND_MAX_SIMTIME"]       = 5e-6
    odict["ENABLE_ORBIT_FOLLOWING"]    = 1
    odict["ENABLE_ORBITWRITE"]         = 1
    odict["ORBITWRITE_MODE"]           = 1
    odict["ORBITWRITE_INTERVAL"]       = 1e-10
    odict["ORBITWRITE_NPOINT"]         = 50002
    odict["ENABLE_MHD"]                = 1

    options.write_hdf5(helpers.testfn, odict, desc="MHD_GCF")


    #**************************************************************************#
    #*           Marker input consisting of an electron and alpha particle     #
    #*                                                                         #
    #**************************************************************************#
    Nmrk   = 2
    ids    = np.array([1, 2])
    weight = np.array([1, 1])
    pitch  = np.array([0.8, 0.9])
    mass   = 1       * np.array([4, m_e_AMU ])
    charge = 1       * np.array([2, -1])
    anum   = 1       * np.array([4, 1])
    znum   = 1       * np.array([1, 1])
    time   = 0       * np.array([1, 1])
    R      = 1       * np.array([0.75, 0.6])
    phi    = 90      * np.array([1, 1])
    z      = 1       * np.array([0, 0])
    zeta   = 2       * np.array([1, 1])
    energy = 10e6    * np.array([.1, 1])
    mrk.write_hdf5(helpers.testfn, Nmrk, ids, mass, charge,
                   R, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc="MHD_GCF")

    #**************************************************************************#
    #*                     Construct a field from geqdsk;                      #
    #*                         add MHD modes                                   #
    #**************************************************************************#


    
    B2DS = eqdsk.generate("input.geq",cocosin=3,desc="MHD_GCF", fn_out=helpers.testfn)
 
    fmhdname = "boozer_maps.out"
    b = Boozermaps(fmhdname)
    b.write_hdf5(helpers.testfn, desc="MHD_GCF")

    #create dummy data
    data = {}
    data["nmode"] = 1 
    data["npsi"] = 500
    data["mmodes"] = np.array([10])
    data["nmodes"] = np.array([11])
    data["amplitude"] = np.array([1e-6])
    data["omega"] = np.array([2e3])
    data["psimin"] = 0
    data["psimax"] = 1
    psi = np.linspace(data["psimin"],data["psimax"],data["npsi"])
    

    data["alpha"] = np.zeros((data["npsi"],data["nmode"]))
    modecenter =0.6
    width = 500.
    for i in range(data["npsi"]):
        data["alpha"][i,0] = np.exp(-width*np.square(psi[i]-modecenter))
    data["phi"] = np.zeros((data["npsi"],data["nmode"]))
    print("Phi fac")
    print(1*data["omega"][0]/c)
    for i in range(data["npsi"]):
        #should be multiplied  by (G+qI)/(nq-m)
        data["phi"][i,0] = (data["omega"][0]/c)*np.exp(-width*np.square(psi[i]-modecenter)) 


    #write that data to file



    mhd.write_hdf5(helpers.testfn, data["nmode"], data["nmodes"],
                   data["mmodes"], data["amplitude"], data["omega"], 
                   data["alpha"], data["phi"],data["npsi"], data["psimin"],
                   data["psimax"], desc="MHD_GCF")


    #**************************************************************************#
    #*                     Rest of the inputs are trivial                      #
    #*                                                                         #
    #**************************************************************************#
    Exyz   = np.array([0, 0, 0])
    E_TC.write_hdf5(helpers.testfn, Exyz, desc="MHD_GCF")

    nwall = 4
    Rwall = np.array([0.1, 100, 100, 0.1])
    zwall = np.array([-100, -100, 100, 100])
    for tname in ["MHD_GCF"]:
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

    This function makes two plots.
    - One that shows conservation of K for all cases
    - And one that shows trajectories on a Rz plane for all cases
    """
    a5 = ascot5.Ascot(helpers.testfn)
    apy = Ascotpy(helpers.testfn)
    apy.init(bfield=True)
    apy.init(boozer=True)
    apy.init(mhd=True) 
    Zmode = np.linspace(-0.70,0.70,500)
    Rmode = np.linspace(0.4,1.1,500)
    nmode = float(a5["MHD_GCF"].mhd.read()["nmodes"][0])
    omega = float(a5["MHD_GCF"].mhd.read()["omega"][0])
    raxis = a5["MHD_GCF"].bfield.read()["axisr"]
    f = plt.figure(figsize=(11.9/2.54, 8/2.54))
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    h1 = f.add_subplot(1,2,1)
    h2 = f.add_subplot(1,2,2)

    colors = ["b", "dodgerblue", "darkgreen", "forestgreen", "r", "tomato"]



    #**************************************************************************#
    #*     Evaluate and plot conservation quantities for MHD_GCF            #
    #*                                                                         #
    #**************************************************************************#
    MHD = {}
    MHD["GCF"] = {}
    orb = a5["MHD_GCF"]["orbit"].read()

    MHD["GCF"]["B"] = np.sqrt(np.power(orb["br"],2) + np.power(orb["bphi"],2) +
                np.power(orb["bz"],2))
    
    MHD["GCF"]["psi"] = apy.evaluate(orb["r"], orb["phi"], orb["z"], orb["time"],"psi")

    # Note that mu is in eV / T

    MHD["GCF"]["time"] = orb["time"]
    MHD["GCF"]["id"]   = orb["id"]
    MHD["GCF"]["r"]    = orb["r"]
    MHD["GCF"]["z"]    = orb["z"]
    MHD["GCF"]["mu"]   = orb["mu"] * e
    MHD["GCF"]["vpar"] = orb["vpar"]
    MHD["GCF"]["charge"] = orb["charge"] 
    
    
    def gamma(mu,mass,vpar, B):
        return np.sqrt( ( 1 + 2 * mu * B / ( mass * c * c) ) /
                        ( 1 - vpar * vpar / ( c * c ) ) )
    def ekin(mu, mass, vpar, B):
        return ( gamma(mu,mass,vpar,B) - 1 )* mass * c * c  

    def ctor(mu, mass, vpar, B, r, charge, psi):
        return gamma(mu,mass,vpar, B) * mass * r * vpar - \
             charge* e*psi   

    def K(mu, mass, vpar, B, r, charge,psi, phi_elec):
        return ekin(mu, mass, vpar, B) + charge *e* phi_elec - \
                    (omega/nmode)*ctor(mu, mass, vpar, B, r, charge,psi)

    MHD["GCF"]["phi"]  = orb["phi"]
    MHD["GCF"]["mhd_phi"] = apy.evaluate(MHD["GCF"]["r"],MHD["GCF"]["phi"],MHD["GCF"]["z"],
                MHD["GCF"]["time"], "mhd_phi")
    id1 = MHD["GCF"]["id"] == 1
    id2 = MHD["GCF"]["id"] == 2



    plot_relerr(h1, MHD["GCF"]["time"][id1], K(MHD["GCF"]["mu"][id1], m_alpha,
                MHD["GCF"]["vpar"][id1], MHD["GCF"]["B"][id1],
                MHD["GCF"]["r"][id1], MHD["GCF"]["charge"][id1], 
                MHD["GCF"]["psi"][id1], MHD["GCF"]["mhd_phi"][id1]), colors[1],
                label="GCF alpha")

    plot_relerr(h1, MHD["GCF"]["time"][id2], K(MHD["GCF"]["mu"][id2], m_e,
                MHD["GCF"]["vpar"][id2], MHD["GCF"]["B"][id2],
                MHD["GCF"]["r"][id2], MHD["GCF"]["charge"][id2], 
                MHD["GCF"]["psi"][id2],MHD["GCF"]["mhd_phi"][id2]), colors[2],
                label="GCF electron")

   
    h2.plot(MHD["GCF"]["r"][id2], MHD["GCF"]["z"][id2],colors[2])
    
    h2.plot(MHD["GCF"]["r"][id1], MHD["GCF"]["z"][id1],colors[1])

  
    apy.plotRz(Rmode,0,Zmode,0,"mhd_br",axes=h2)

  #**************************************************************************#
    #*                 Finalize and print and show the figure                  #
    #*                                                                         #
    #**************************************************************************#

    h1.set_xlim(0, 5e-6)
    h1.xaxis.set(ticks=[0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6], ticklabels=[0,1,2,3,4,5])
    h1.tick_params(axis='y', direction='out')
    h1.tick_params(axis='x', direction='out')
    h1.spines['right'].set_visible(False)
    h1.spines['top'].set_visible(False)
    h1.yaxis.set_ticks_position('left')
    h1.xaxis.set_ticks_position('bottom')
    h1.set(ylabel="$\Delta K/K_0$",xlabel="Time[10-6 s]")
    h1.legend()

    h2.set(xlabel="$R$ [m]", ylabel = "$z$ [m]")

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
