import numpy as np
import subprocess
import sys
from os import chdir 
from os.path import dirname, realpath, sep, pardir 
sys.path.append(dirname(realpath(__file__)) + sep + pardir + sep + 'ui')
import scipy.constants as constants
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.options as options
import a5py.ascot5io.B_GS as B_GS
import a5py.ascot5io.markers as markers
import analyticBKGpsifun as psifun
import testfunctions as tf
from testcase import createbase

def run():

    ELEMENTARY_CHARGE = constants.elementary_charge
    AMU2KG = constants.physical_constants['atomic mass unit-kilogram relationship'][0]
    BOLTZMANN_CONSTANT = constants.Boltzmann
    EV2K = constants.physical_constants['electron volt-kelvin relationship'][0]

    # Test options
    simmode = [2, 2, 1]
    adaptive = [0, 1, 0]
    timestep = [1e-8, 1e-8, 1.e-10]
    tolorb = 1e-8
    tolcol = 1e-2
    writedt = 2e-7
    simtime = 1e-3
    Nmrk = 4

    # Proton
    m = 1.00727647
    q = 1

    # Init testbase.
    fn = ['GCfixed.h5', 'GCadaptive.h5', 'GO.h5'] #@GO adaptive == fixed
    Bxyz = np.array([1, 0, 0])
    Exyz = np.array([0, 0, 0])
    n = 1e20
    T = 1e3

    # Init magnetic field.
    R0 = 6.2
    z0 = 0
    B_phi0 = 5.3
    psi_mult = 200
    psi_coeff = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 
                          3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 
                          7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, 
                          -0.155])

    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = m*np.ones(ids.shape)
    charge = q*np.ones(ids.shape)
    R      = np.linspace(6.4,7.4,Nmrk)
    phi    = 0*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    energy = 1.0e5*np.ones(ids.shape)
    pitch  = np.array([-0.9, -0.3, 0.3, 0.9])
    theta  = 0*np.ones(ids.shape)

    for i in range(0,len(fn)):
        createbase(fn[i], Bxyz, Exyz, n, T)
        B_GS.write_hdf5(fn[i], R0, z0, B_phi0, psi_mult, psi_coeff)

        # Options
        o = options.read_hdf5(fn[i])
        o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode[i]
        o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive[i]
        o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
        o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep[i]
        o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
        o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
        o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
        o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
        o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
        o["ENDCOND_SIMTIMELIM"]        = 0*o["ENDCOND_SIMTIMELIM"] + 1
        o["ENDCOND_MAX_SIM_TIME"]      = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime
        o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + writedt
        options.write_hdf5(fn[i],o)

        markers.write_hdf5_guidingcenters(fn[i], Nmrk, ids, mass, charge, R, phi, z, energy, pitch, theta, weight, time)

        # Simulate.
        subprocess.call(["./ascot5_main", "--in="+fn[i][0:-3]])

    # Read orbits
    orb, t, B, R, z, e, mu, pphi = ([] for i in range(8))
    for i in range(0,len(fn)):
        if simmode[i] == 1:
            orb.append(ascot5.read_hdf5(fn[i],"orbits")["orbits"]["fo"])
        else:
            orb.append(ascot5.read_hdf5(fn[i],"orbits")["orbits"]["gc"])

        t.append(orb[i]["time"])
        B.append(np.sqrt(np.power(orb[i]["B_R"],2) + np.power(orb[i]["B_phi"],2) + np.power(orb[i]["B_z"],2)))
        R.append(orb[i]["R"])
        z.append(orb[i]["z"])
        if simmode[i] == 1:
            v = np.sqrt(np.power(orb[i]["v_R"],2) + np.power(orb[i]["v_phi"],2) + np.power(orb[i]["v_z"],2))
            e.append(0.5 * m * np.power(v,2) * AMU2KG / ELEMENTARY_CHARGE)
            
            pitch = (orb[i]["B_R"] * orb[i]['v_R'] + orb[i]["B_phi"] * orb[i]['v_phi'] + orb[i]["B_z"] * orb[i]['v_z']) / ( v * B[i] )

            mu.append(0.5*(1-np.power(pitch,2))*np.power(v,2)*m*AMU2KG/(B[i]*ELEMENTARY_CHARGE))
        else:
            e.append(orb[i]["mu"] * B[i] + 0.5 * m * np.power(orb[i]["vpar"],2) * AMU2KG / ELEMENTARY_CHARGE)
            mu.append(orb[i]["mu"])
    
        #Canonical momentum
        c = psi_coeff
        psi = psi_mult*psifun.psi0(R[i]/R0,z[i]/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])
        if simmode[i] == 2:
            vtor = orb[i]["vpar"] * orb[i]["B_phi"] / B[i]
            pphi.append(mass[0]* AMU2KG * R[i] * vtor + charge[0]* ELEMENTARY_CHARGE*psi)
                #tf.canonicalMomentum(mass[0]* AMU2KG,R[i],orb[i]["vpar"],charge[0]* ELEMENTARY_CHARGE,psi[i]))
        else:
            pphi.append(mass[0]* AMU2KG * R[i] * orb[i]["v_phi"] + charge[0]* ELEMENTARY_CHARGE*psi)
                #tf.canonicalMomentum(mass[0]* AMU2KG,R[i],vpar,charge[0]* ELEMENTARY_CHARGE,psi[i]))

    # Plot if needed.
    
    
    plt.figure()
    for i in range(len(fn)):
        plt.plot(t[i],e[i],'.')
    plt.show()

    plt.figure()
    for i in range(len(fn)):
        plt.plot(t[i],mu[i],'.')
    plt.show()

    plt.figure()
    for i in range(len(fn)):
        plt.plot(R[i],z[i],'.')
    plt.show()
    
    plt.figure()
    for i in range(len(fn)):
        plt.plot(t[i],pphi[i],'.')
    plt.show()
    

    # Write relevant data in files for plotting in Matlab

    #canonical momentum
    with open('pphi.txt', 'w') as f:
        for run in pphi:
            for val in run:
                f.write("%s," % val)
            f.write("\n")

    #time
    with open('time.txt', 'w') as f:
        for run in t:
            for val in run:
                f.write("%s," % val)
            f.write("\n")

    #energy
    with open('energy.txt', 'w') as f:
        for run in e:
            for val in run:
                f.write("%s," % val)
            f.write("\n")

    #magnetic moment
    with open('mu.txt', 'w') as f:
        for run in mu:
            for val in run:
                f.write("%s," % val)
            f.write("\n")

    #location in R-coordinate
    with open('R.txt', 'w') as f:
        for run in R:
            for val in run:
                f.write("%s," % val)
            f.write("\n")

    #location in z-coordinate
    with open('z.txt', 'w') as f:
        for run in z:
            for val in run:
                f.write("%s," % val)
            f.write("\n")

    # Compare.


    # Clean.
    #for f in fn:
    #    subprocess.call(["rm", f])

if __name__ == '__main__':
    run()
