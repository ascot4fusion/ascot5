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
import a5py.postprocessing.physicslib as physlib
import a5py.preprocessing.analyticequilibrium as psifun
import testfunctions as tf
from testcase import createbase

ELEMENTARY_CHARGE = constants.elementary_charge
AMU2KG = constants.physical_constants['atomic mass unit-kilogram relationship'][0]
BOLTZMANN_CONSTANT = constants.Boltzmann
EV2K = constants.physical_constants['electron volt-kelvin relationship'][0]

psi_coeff = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 
                      3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 
                      7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, 
                      -0.155])
psi_mult = 200
R0 = 6.2
z0 = 0

simmode = [2, 2, 1]
fn = ['orbfol_GCFIX.h5', 'orbfol_GCADA.h5', 'orbfol_GO.h5']

def init(fast):
    

    # Test options
   
    adaptive      = [0, 1, 0]
    timestep      = [1e-11, 1e-11, 1e-13]
    tolorb        = 1e-9
    Nmrk          = 4

    if(fast):
        simtime = 5e-6
    else:
        simtime = 1e-4
    writedt       = simtime/1e5

    # Proton and electron info (mass [amu], charge [e], energy [eV])
    mp = 1.00727647
    qp = 1
    Ep = 1e6

    me = 5.485799090e-4
    qe = -1
    Ee = 1e8


    # Init testbase.
    Bxyz = np.array([1, 0, 0])
    Exyz = np.array([0, 0, 0])
    n = 1e20
    T = 1e3

    # Init magnetic field.
    B_phi0 = 5.3
    
    
    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = me*np.ones(ids.shape)
    charge = qe*np.ones(ids.shape)
    R      = 7.6*np.ones(ids.shape)
    phi    = 0*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    energy = 1e7*np.ones(ids.shape)
    pitch  = np.array([-0.9, -0.5, 0.9, 0.5])
    theta  = 0*np.ones(ids.shape)

    for i in range(0,len(fn)):
        createbase(fn[i], Bxyz, Exyz, n, T)
        B_GS.write_hdf5(fn[i], R0, z0, B_phi0, psi_mult, psi_coeff)

        # Options
        o = options.read_hdf5(fn[i])
        o["SIM_MODE"]                      = 0*o["SIM_MODE"] + simmode[i]
        o["ENABLE_ADAPTIVE"]               = 0*o["ENABLE_ADAPTIVE"] + adaptive[i]
        o["FIXEDSTEP_USE_USERDEFINED"]     = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
        o["FIXEDSTEP_USERDEFINED"]         = 0*o["FIXEDSTEP_USERDEFINED"] + timestep[i]
        o["ADAPTIVE_TOL_ORBIT"]            = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
        o["ENABLE_ORBIT_FOLLOWING"]        = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
        o["ENABLE_ORBITWRITE"]             = 0*o["ENABLE_ORBITWRITE"] + 1
        o["ORBITWRITE_MODE"]               = 0*o["ORBITWRITE_MODE"] + 1
        o["ENDCOND_SIMTIMELIM"]            = 0*o["ENDCOND_SIMTIMELIM"] + 1
        o["ENDCOND_MAX_SIM_TIME"]          = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime
        o["ORBITWRITE_INTERVAL"]           = 0*o["ORBITWRITE_INTERVAL"] + writedt
        options.write_hdf5(fn[i],o)

        markers.write_hdf5_guidingcenters(fn[i], Nmrk, ids, mass, charge, R, phi, z, energy, pitch, theta, weight, time)

def run():
    # Simulate.
    for i in range(len(fn)):
        subprocess.call(["./ascot5_main", "--in="+fn[i][0:-3]])

def check(plot):
    Etol = np.array([1e-6, 1e-6, 1e-6])
    mutol = np.array([1e-8, 1e-8, 5e-2])
    cpphitol = np.array([1e-3, 1e-3, 1e-3])

    # Read orbits
    ids, t, R, z, Ekin, mu, cpphi = ([] for i in range(7))
    for i in range(0,len(fn)):
        if simmode[i] == 1:
            orb = ascot5.read_hdf5(fn[i],"orbits")["orbits"]["fo"]
        else:
            orb = ascot5.read_hdf5(fn[i],"orbits")["orbits"]["gc"]

        ids.append(orb["id"])
        t.append(orb["time"])
        R.append(orb["R"])
        z.append(orb["z"])

        c = psi_coeff
        psi = psi_mult*psifun.psi0(R[i]/R0,z[i]/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])
        B = np.sqrt(np.power(orb["B_R"],2) + np.power(orb["B_phi"],2) + np.power(orb["B_z"],2))

        if simmode[i] == 1:
            v = np.sqrt(np.power(orb["v_R"],2) + np.power(orb["v_phi"],2) + np.power(orb["v_z"],2))
            pitch = (orb["B_R"] * orb['v_R'] + orb["B_phi"] * orb['v_phi'] + orb["B_z"] * orb['v_z']) / (v*B)
            gamma = 1/np.sqrt( 1-np.power(v/constants.c,2) )

            Ekin.append( (gamma-1)*orb["mass"]*AMU2KG * constants.c * constants.c / ELEMENTARY_CHARGE)
            mu.append(0.5*(1-np.power(pitch,2))*np.power(v*gamma,2)*orb["mass"]*AMU2KG/(B*ELEMENTARY_CHARGE))
            cpphi.append(gamma * orb["mass"] * AMU2KG * R[i] * orb["v_phi"] + orb["charge"] * ELEMENTARY_CHARGE * psi)
        else:
            Ekin.append(physlib.Ekin(massamu=orb["mass"], mueVperT=orb["mu"], vpar=orb["vpar"], Btot=B)/ELEMENTARY_CHARGE)
            gamma = 1 + Ekin[i]/(orb["mass"]*AMU2KG * constants.c * constants.c / ELEMENTARY_CHARGE)

            mu.append(orb["mu"])
            vtor = orb["vpar"] * orb["B_phi"] / B
            cpphi.append(gamma * orb["mass"] * AMU2KG * R[i] * vtor + orb["charge"] * ELEMENTARY_CHARGE * psi)
            

    # Compare.
    isOkay = np.ones((3,4),dtype='int32')
    Ei = np.zeros(isOkay.shape)
    Ef = np.zeros(isOkay.shape)
    Ee = np.zeros(isOkay.shape)
    mui = np.zeros(isOkay.shape)
    muf = np.zeros(isOkay.shape)
    mue = np.zeros(isOkay.shape)
    cpphii = np.zeros(isOkay.shape)
    cpphif = np.zeros(isOkay.shape)
    cpphie = np.zeros(isOkay.shape)
    
    for i in range(0,len(fn)):
        for j in range(0,4):
            idx = ids[i] == j+1
            E = Ekin[i][idx]
            Ei[i][j] = np.max(E)
            Ef[i][j] = np.min(E)
            Ee[i][j] = np.absolute((Ei[i][j]-Ef[i][j])/Ei[i][j])

            muu = mu[i][idx]
            mui[i][j] = np.max(muu)
            muf[i][j] = np.min(muu)
            mue[i][j] = np.absolute((mui[i][j]-muf[i][j])/mui[i][j])

            cpphix = cpphi[i][idx]
            cpphii[i][j] = np.max(cpphix)
            cpphif[i][j] = np.min(cpphix)
            cpphie[i][j] = np.absolute((cpphii[i][j]-cpphif[i][j])/cpphii[i][j])
            
            if( Ee[i][j]  > Etol[i]):
                isOkay[i,j] = 0
            if( mue[i][j] > mutol[i]):
                isOkay[i,j] = 0
            if( cpphie[i][j]  > cpphitol[i]):
                isOkay[i,j] = 0

    ok = ["WRONG!", "OK"]

    print("\nInitial and final energy, mag. mom., and can. tor. mom. and error\n")
    mode = ["GCA", "GCF", "GO"]
    for j in range(2,-1,-1):
        for i in range(0,4):
            print("\n" + mode[j] + " " + str(i+1) +
                  " E: " + str(Ei[j][i]) + ", " + str(Ef[j][i]) + ", " + str(Ee[j][i]) + 
                  "\n mu: " + str(mui[j][i]) + ", " + str(muf[j][i]) + ", " + str(mue[j][i]) + 
                  "\n cpphi: " + str(cpphii[j][i]) + ", " + str(cpphif[j][i]) + ", " + str(cpphie[j][i]) + 
                  "\n " + ok[isOkay[j][i]])

    isOkay = np.sum(np.sum(isOkay))
    if(isOkay == 3*4):
        print("Test succeeded.")
    else:
        print("Test failed!")

    # Plot if needed.
    if(plot):
        c = ['b', 'g', 'r']
        f  = plt.figure()
        a1 = f.add_subplot(2,2,1)
        a2 = f.add_subplot(2,2,2)
        a3 = f.add_subplot(2,2,3)
        a4 = f.add_subplot(2,2,4)

        for i in range(2,-1,-1):
            for j in range(0,4):
                idx = ids[i] == j+1
                Ei = Ekin[i][idx]
                ti = t[i][idx]
                mui = mu[i][idx]
                cpphii = cpphi[i][idx]
                Ri = R[i][idx]
                zi = z[i][idx]

                a1.plot(ti,Ei/Ei[0],color=c[i])
                a2.plot(ti,mui/mui[0],color=c[i])
                a3.plot(ti,cpphii/cpphii[0],color=c[i])
                a4.plot(Ri,zi,color=c[i])

        a1.set_xlabel("Time (s)")
        a1.set_ylabel("Relative change in energy")

        a2.set_xlabel("Time (s)")
        a2.set_ylabel("Relative change in mu")

        a3.set_xlabel("Time (s)")
        a3.set_ylabel("Relative change in cpphhi")

        a4.set_xlabel("R (m)")
        a4.set_ylabel("z (m)")
        
        plt.show()

    return isOkay

def clean():
    # Clean
    for f in fn:
        subprocess.call(["rm", f])

if __name__ == '__main__':
    print("\n\
    This test validates the orbit-following.\n\
    Markers are simulated in ITER-like background without collisions.\n\
    Markers are relativistic electrons (Ekin = 100 MeV). Simulations \n\
    are repeated with GO, GC fixed and GC adaptive modes.\n\
    \n\
    In GC modes, energy and canonical toroidal momentum may drift but this\n\
    drift should decrease with time step. Magnetic moment should be conserved\n\
    exactly.\n\
    \n\
    In GO mode, energy is conserved exactly but canonical momentum and\n\
    magnetic moment might drift.\n\
    ")

    init(True)
    run()
    check(True)
    clean()
