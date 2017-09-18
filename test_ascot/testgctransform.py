import numpy as np
import subprocess
import sys
import scipy.constants as constants
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.options as options
import a5py.ascot5io.B_GS as B_GS
import a5py.ascot5io.markers as markers
import a5py.postprocessing.physicslib as physlib
import testfunctions as tf
from testcase import createbase
   

def run():
    ELEMENTARY_CHARGE = constants.elementary_charge
    AMU2KG = constants.physical_constants['atomic mass unit-kilogram relationship'][0]
    BOLTZMANN_CONSTANT = constants.Boltzmann
    EV2K = constants.physical_constants['electron volt-kelvin relationship'][0]

    # Test options
    writedt = 1e-11
    simtime = 4e-6

    # Proton
    m = 1.00727647
    q = 1

    # Init testbase.
    fn = ['GO.h5', 'GO2GC.h5', 'GC.h5']
    Bxyz = np.array([0, 1, 0])
    Exyz = np.array([0, 0, 0])
    n = 1e20
    T = 1e3
    createbase(fn[0], Bxyz, Exyz, n, T)

    # Init magnetic field.
    R0 = 6.2
    z0 = 0
    B_phi0 = 5.3
    psi_mult = 200
    psi_coeff = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 
                          3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 
                          7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, 
                          -0.155])
    B_GS.write_hdf5(fn[0], R0, z0, B_phi0, psi_mult, psi_coeff)

    # Markers
    Nmrk   = 1
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = m*np.ones(ids.shape)
    charge = q*np.ones(ids.shape)
    R      = 7.4*np.ones(ids.shape)
    phi    = 0*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    energy = 1.0e3*np.ones(ids.shape)
    pitch  = 0.999*np.ones(ids.shape)
    theta  = 0*np.ones(ids.shape)  

    # Options
    o = options.read_hdf5(fn[0])
    o["SIM_MODE"]                  = 0*o["SIM_MODE"] + 1
    o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
    o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + 1e-11
    o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
    o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
    o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
    o["ENDCOND_SIMTIMELIM"]        = 0*o["ENDCOND_SIMTIMELIM"] + 1
    o["ENDCOND_MAX_SIM_TIME"]      = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime
    o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + writedt
    options.write_hdf5(fn[0],o)
    
    markers.write_hdf5_guidingcenters(fn[0], Nmrk, ids, mass, charge, R, phi, z, energy, pitch, theta, weight, time)
    
    subprocess.call(["cp", fn[0], fn[1]])
    subprocess.call(["cp", fn[0], fn[2]])

    o = options.read_hdf5(fn[1])
    o["RECORD_GO_AS_GC"] = 0*o["RECORD_GO_AS_GC"] + 1
    options.write_hdf5(fn[1],o)

    o = options.read_hdf5(fn[2])
    o["SIM_MODE"]                  = 0*o["SIM_MODE"] + 2
    o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + 1e-9
    options.write_hdf5(fn[2],o)

    # Simulate.
    subprocess.call(["./ascot5_main", "--in="+fn[0][0:-3]])
    subprocess.call(["./ascot5_main", "--in="+fn[1][0:-3]])
    subprocess.call(["./ascot5_main", "--in="+fn[2][0:-3]])

    
    go    = ascot5.read_hdf5(fn[0],"orbits")["orbits"]["fo"]
    go2gc = ascot5.read_hdf5(fn[1],"orbits")["orbits"]["gc"]
    gc    = ascot5.read_hdf5(fn[2],"orbits")["orbits"]["gc"]

    pitch = physlib.pitch(massamu=go["mass"], vR=go["v_R"], vphi=go["v_phi"], vz=go["v_z"],
                          BR=go["B_R"], Bphi=go["B_phi"], Bz=go["B_z"])
    gamma = physlib.gamma(massamu=go["mass"], vR=go["v_R"], vphi=go["v_phi"], vz=go["v_z"],
                          BR=go["B_R"], Bphi=go["B_phi"], Bz=go["B_z"])
    govpar = pitch * physlib.vecnorm(go["v_R"],go["v_phi"],go["v_z"])
    gomu   = gamma*gamma*(1-pitch*pitch) * go["mass"] * AMU2KG * np.power(physlib.vecnorm(go["v_R"],go["v_phi"],go["v_z"]),2)/physlib.vecnorm(go["B_R"],go["B_phi"],go["B_z"])
    gomu   = 0.5*gomu/ELEMENTARY_CHARGE

    if True:
        plt.figure()
        plt.plot(go["R"], go["z"])
        plt.plot(go2gc["R"], go2gc["z"])
        plt.plot(gc["R"], gc["z"])
    else:
        plt.figure()
        x,y = pol2cart(go["R"], go["phi"]*np.pi/180)
        plt.plot(x,go["z"])
        x,y = pol2cart(go2gc["R"], go2gc["phi"]*np.pi/180)
        plt.plot(x,go2gc["z"])
        x,y = pol2cart(gc["R"], gc["phi"]*np.pi/180)
        plt.plot(x,gc["z"])
        plt.axis('equal')
        


    plt.figure()
    plt.plot(go["time"], gomu)
    plt.plot(go2gc["time"], go2gc["mu"])
    plt.plot(gc["time"], gc["mu"])

    plt.figure()
    plt.plot(go["time"], govpar)
    plt.plot(go2gc["time"], go2gc["vpar"])
    plt.plot(gc["time"], gc["vpar"])

    plt.show()

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

if __name__ == '__main__':
    run()
