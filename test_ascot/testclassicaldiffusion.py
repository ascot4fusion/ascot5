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
    simmode = 1
    adaptive = 1
    timestep = 1.e-10
    tolorb = 1e-8
    tolcol = 1e-2
    orbfollowing = 1
    simtime = 1e-4
    Nmrk = 100

    # Proton
    m = 1.00727647
    q = 1

    # Init testbase.
    fn = "test1.h5"
    Bxyz = np.array([0, 1, 0])
    Exyz = np.array([0, 0, 0])
    n = 1e20
    T = 1e4
    createbase(fn, Bxyz, Exyz, n, T)

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
    o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 0
    o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
    o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
    o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + 1e-7
    options.write_hdf5(fn,o)

    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = m*np.ones(ids.shape)
    charge = q*np.ones(ids.shape)
    R      = 5*np.ones(ids.shape)
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
        x = orb["R"] * np.cos(np.deg2rad(orb["phi"]))
        y = orb["R"] * np.sin(np.deg2rad(orb["phi"]))
        z = orb["z"]
        t = orb["time"]

        
        
    else:
        orb = ascot5.read_hdf5(fn,"orbits")["orbits"]["fo"]
        x = orb["R"] * np.cos(np.deg2rad(orb["phi"]))
        y = orb["R"] * np.sin(np.deg2rad(orb["phi"]))
        z = orb["z"]
        t = orb["time"]
    
    r = np.sqrt(np.power(x,2) + np.power(z,2))

    state = ascot5.read_hdf5(fn,"states")["states"]["endstate"]
    xf = state["R"] * np.cos(np.deg2rad(state["phi"]))
    yf = state["R"] * np.sin(np.deg2rad(state["phi"]))
    zf = state["z"]
    rf = np.sqrt(np.power(xf,2) + np.power(zf,2))
    muf = state["mu"]

    state = ascot5.read_hdf5(fn,"states")["states"]["inistate"]
    xi = state["R"] * np.cos(np.deg2rad(state["phi"]))
    yi = state["R"] * np.sin(np.deg2rad(state["phi"]))
    zi = state["z"]
    ri = np.sqrt(np.power(xi,2) + np.power(zi,2))
    mui = state["mu"]

    # Plot if needed.
    #plt.plot(t,r,'.')
    #plt.show()

    plt.figure()
    n, bins, patches = plt.hist(muf-mui, 50, normed=1, facecolor='green', alpha=0.75)
    plt.show()


    # Compare.
    #D = np.mean(np.power(ri-rf,2) / (2*simtime))
    D = np.var(ri-rf)/simtime
    print("D = " + ("%.5f" % D))

    # Clean.
    subprocess.call(["rm", fn])

if __name__ == '__main__':
    run()

