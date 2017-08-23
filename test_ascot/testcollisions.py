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
    timestep = 1.e-6
    tolorb = 1e-8
    tolcol = 1e-2
    orbfollowing = 0
    Nmrk = 100

    # Proton
    m = 1.00727647
    q = 1

    # Init testbase.
    fn = "test1.h5"
    Bxyz = np.array([1, 0, 0])
    Exyz = np.array([0, 0, 0])
    n = 1e20
    T = 1e3
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

    # Common options
    o = options.read_hdf5(fn)
    o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode
    o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive
    o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
    o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep
    o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
    o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
    o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + orbfollowing
    o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 1
    o["ENABLE_RZVparaVperp_DIST"]  = 0*o["ENABLE_RZVparaVperp_DIST"] + 1
    o["DIST_RZVparaVperp_MIN_R"]   = 0*o["DIST_RZVparaVperp_MIN_R"] + 4
    o["DIST_RZVparaVperp_MAX_R"]   = 0*o["DIST_RZVparaVperp_MAX_R"] + 8
    o["DIST_RZVparaVperp_BIN_R"]   = 0*o["DIST_RZVparaVperp_BIN_R"] + 1
    o["DIST_RZVparaVperp_MIN_Z"]   = 0*o["DIST_RZVparaVperp_MIN_Z"] - 5
    o["DIST_RZVparaVperp_MAX_Z"]   = 0*o["DIST_RZVparaVperp_MAX_Z"] + 5
    o["DIST_RZVparaVperp_BIN_Z"]   = 0*o["DIST_RZVparaVperp_BIN_Z"] + 1
    o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
    o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
    options.write_hdf5(fn,o)

    # Clone base into two copies.
    fn1 = fn
    fn2 = "test2.h5"
    subprocess.call(["cp",fn,fn2])

    # Slowing down options and markers
    o = options.read_hdf5(fn1)
    o["ENDCOND_ENERGYLIM"]                = 0*o["ENDCOND_ENERGYLIM"] + 1
    o["ENDCOND_MIN_ENERGY"]               = 0*o["ENDCOND_MIN_ENERGY"] + 1e3
    o["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = 0*o["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] + 0.5
    o["ORBITWRITE_INTERVAL"]              = 0*o["ORBITWRITE_INTERVAL"] + 1e-4
    o["DIST_RZVparaVperp_MIN_VPARA"]      = 0*o["DIST_RZVparaVperp_MIN_VPARA"] - 7e7
    o["DIST_RZVparaVperp_MAX_VPARA"]      = 0*o["DIST_RZVparaVperp_MAX_VPARA"] + 7e7
    o["DIST_RZVparaVperp_BIN_VPARA"]      = 0*o["DIST_RZVparaVperp_BIN_VPARA"] + 100
    o["DIST_RZVparaVperp_MIN_VPERP"]      = 0*o["DIST_RZVparaVperp_MIN_VPERP"] + 0
    o["DIST_RZVparaVperp_MAX_VPERP"]      = 0*o["DIST_RZVparaVperp_MAX_VPERP"] + 7e7
    o["DIST_RZVparaVperp_BIN_VPERP"]      = 0*o["DIST_RZVparaVperp_BIN_VPERP"] + 50
    options.write_hdf5(fn1,o)

    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = m*np.ones(ids.shape)
    charge = q*np.ones(ids.shape)
    R      = 6.4*np.ones(ids.shape)
    phi    = 0*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    theta  = 0*np.ones(ids.shape)
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    
    energy = 1.0e5*np.ones(ids.shape)
    pitch  = 0.9*np.ones(ids.shape)
    markers.write_hdf5_guidingcenters(fn1, Nmrk, ids, mass, charge, 
                                      R, phi, z, energy, pitch, theta, weight, time)

    # Thermal options and markers
    o = options.read_hdf5(fn2)
    o["ENDCOND_SIMTIMELIM"]          = 0*o["ENDCOND_SIMTIMELIM"] + 1
    o["ENDCOND_MAX_SIM_TIME"]        = 0*o["ENDCOND_MAX_SIM_TIME"] + 1e-4
    o["ORBITWRITE_INTERVAL"]         = 0*o["ORBITWRITE_INTERVAL"] + 1e-5
    o["DIST_RZVparaVperp_MIN_VPARA"] = 0*o["DIST_RZVparaVperp_MIN_VPARA"] - 1e7
    o["DIST_RZVparaVperp_MAX_VPARA"] = 0*o["DIST_RZVparaVperp_MAX_VPARA"] + 1e7
    o["DIST_RZVparaVperp_BIN_VPARA"] = 0*o["DIST_RZVparaVperp_BIN_VPARA"] + 100
    o["DIST_RZVparaVperp_MIN_VPERP"] = 0*o["DIST_RZVparaVperp_MIN_VPERP"] + 0
    o["DIST_RZVparaVperp_MAX_VPERP"] = 0*o["DIST_RZVparaVperp_MAX_VPERP"] + 1e7
    o["DIST_RZVparaVperp_BIN_VPERP"] = 0*o["DIST_RZVparaVperp_BIN_VPERP"] + 50
    options.write_hdf5(fn2,o)

    energy = 1.0e3*np.ones(ids.shape)
    pitch  = 0.9*np.ones(ids.shape)
    markers.write_hdf5_guidingcenters(fn2, Nmrk, ids, mass, charge, 
                                      R, phi, z, energy, pitch, theta, weight, time)


    # Simulate.
    subprocess.call(["./ascot5_main", "--in="+fn1[0:-3]])
    subprocess.call(["./ascot5_main", "--in="+fn2[0:-3]])

    if simmode == 2:
        # Read distribution and orbits from slowing down
        is1 = ascot5.read_hdf5(fn1,"orbits")["orbits"]["gc"]
        t1 = is1["time"]
        B  = np.sqrt(np.power(is1["B_R"],2) + np.power(is1["B_phi"],2) + np.power(is1["B_z"],2))
        e1 = is1["mu"] * B + 0.5*m*np.power(is1["vpar"],2) * AMU2KG / ELEMENTARY_CHARGE
        
        # Read endstate and orbits from thermal
        is2 = ascot5.read_hdf5(fn2,"orbits")["orbits"]["gc"]
        t2 = is2["time"]
        B  = np.sqrt(np.power(is2["B_R"],2) + np.power(is2["B_phi"],2) + np.power(is2["B_z"],2))
        e2 = is2["mu"] * B + 0.5*m*np.power(is2["vpar"],2) * AMU2KG / ELEMENTARY_CHARGE

        is2 = ascot5.read_hdf5(fn2,"states")["states"]["endstate"]
        B  = np.sqrt(np.power(is2["B_R"],2) + np.power(is2["B_phi"],2) + np.power(is2["B_z"],2))
        thermaldist = is2["mu"] * B + 0.5*m*np.power(is2["vpar"],2) * AMU2KG / ELEMENTARY_CHARGE
    else:
        # Read distribution and orbits from slowing down
        is1 = ascot5.read_hdf5(fn1,"orbits")["orbits"]["fo"]
        t1 = is1["time"]
        e1 = 0.5 * m * ( np.power(is1["v_R"],2) + np.power(is1["v_phi"],2) + np.power(is1["v_z"],2) )* AMU2KG / ELEMENTARY_CHARGE
        
        # Read endstate and orbits from thermal
        is2 = ascot5.read_hdf5(fn2,"orbits")["orbits"]["fo"]
        t2 = is2["time"]
        e2 = 0.5 * m * ( np.power(is2["v_R"],2) + np.power(is2["v_phi"],2) + np.power(is2["v_z"],2) )* AMU2KG / ELEMENTARY_CHARGE
    
    # Plot if needed.
    
    """
    plt.figure()
    plt.plot(t1,e1)
    plt.show()

    plt.figure()
    plt.plot(t2,e2)
    plt.show()
    """
    E_edges = np.linspace(0,2,1000)*1e5
    xi_edges = np.linspace(-1,1,20)
    RzExi = disttrans.vpavpe2Epitch(ascot5.read_hdf5(fn1,"dists")["dists"], E_edges, xi_edges, m)
    RzExi = np.squeeze(RzExi['ordinate'])

    plt.figure()
    Edist = RzExi.sum(axis=1)
    plt.plot(Edist)
    plt.show()

    plt.figure()
    n, bins, patches = plt.hist(thermaldist, 50, normed=1, facecolor='green', alpha=0.75)
    plt.show()
    # Compare.


    # Clean.
    subprocess.call(["rm", fn1])
    subprocess.call(["rm", fn2])

if __name__ == '__main__':
    run()
