import numpy as np
import subprocess
import scipy.constants as constants
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.options as options
import a5py.ascot5io.B_GS as B_GS
import a5py.ascot5io.markers as markers
import a5py.postprocessing.disttrasnformations as disttrans
from testcase import createbase

def run():

    ELEMENTARY_CHARGE = constants.elementary_charge
    AMU2KG = 1.66053904e-27

    # Test options
    simmode = 2
    adaptive = 1
    timestep = 1.e-8
    tolorb = 1e-8
    tolcol = 1e-2
    orbfollowing = 1
    simtime = 2e-2
    Nmrk = 1000

    # Proton
    m = 1.00727647
    q = 1

    # Init testbase.
    fn = "test1.h5"
    Bxyz = np.array([0, 1, 0])
    Exyz = np.array([0.0, 0.0, 0.0])
    n = 1e20
    T = 1e4
    createbase(fn, Bxyz, Exyz, n, T)

    # Init magnetic field.
    R0 = 6.2
    z0 = 0
    B_phi0 = 5.3
    psi_mult = 200
    psi_coeff = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 
                          3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 
                          7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, 
                          -0.155])

    B_GS.write_hdf5(fn, R0, z0, B_phi0, psi_mult, psi_coeff)


    # Options
    o = options.read_hdf5(fn)
    o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode
    o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive
    o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
    o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep
    o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
    o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
    o["ENDCOND_SIMTIMELIM"]        = 0*o["ENDCOND_SIMTIMELIM"] + 1
    o["ENDCOND_MAX_SIM_TIME"]      = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime
    o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + orbfollowing
    o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 1
    o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
    o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
    o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + 1e-6
    options.write_hdf5(fn,o)

    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = m*np.ones(ids.shape)
    charge = q*np.ones(ids.shape)
    R      = 7.8*np.ones(ids.shape)
    phi    = 0*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    theta  = np.random.rand(1,Nmrk)*2*np.pi
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    
    energy = 1.0e3*np.ones(ids.shape)
    pitch  = np.linspace(-1,1,Nmrk)
    markers.write_hdf5_guidingcenters(fn, Nmrk, ids, mass, charge, 
                                      R, phi, z, energy, pitch, theta, weight, time)


    # Simulate.
    subprocess.call(["./ascot5_momcoll", "--in="+fn[0:-3]])

    if simmode == 2:
        orb = ascot5.read_hdf5(fn,"orbits")["orbits"]["gc"]   
        
    else:
        orb = ascot5.read_hdf5(fn,"orbits")["orbits"]["fo"]
        
    t   = orb["time"]
    rho = orb["rho"]

    state = ascot5.read_hdf5(fn,"states")["states"]["endstate"]
    rhof = state["rho"]

    state = ascot5.read_hdf5(fn,"states")["states"]["inistate"]
    rhoi = state["rho"]

    # Plot if needed.
    #plt.plot(t,r,'.')
    #plt.show()

    plt.figure()
    n, bins, patches = plt.hist(rhof, 50, normed=1, facecolor='green', alpha=0.75)
    plt.show()


    # Compare.
    #D = np.mean(np.power(ri-rf,2) / (2*simtime))
    D = np.var(rhoi-rhof)/simtime
    print("D = " + ("%.5f" % D))

    # Clean.
    subprocess.call(["rm", fn])

if __name__ == '__main__':
    run()

