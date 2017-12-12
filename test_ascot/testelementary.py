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

fn = ["elemGO.h5", "elemGCF.h5", "elemGCA.h5"]

ELEMENTARY_CHARGE = constants.elementary_charge
AMU2KG = 1.66053904e-27

simtime = 2e-5

mprt = 1.00727647*2
qprt = 1*2

Temp = 1e3
dens = 1e22
Bnorm = 5
Eprt = 3.5e6

simmode = [1, 2, 2]

def init():

    

    # Test options
    adaptive = [0, 0, 1]
    timestep = [5.e-9, 1e-8, 1e-8]
    tolorb = 1e-8
    tolcol = 1e-1
    Nmrk = 1



    # Init testbase.
    Bxyz = np.array([0, Bnorm, 0])
    Exyz = np.array([0, 0, 0])
    createbase(fn[0], Bxyz, Exyz, dens, Temp)

    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = mprt*np.ones(ids.shape)
    charge = -qprt*np.ones(ids.shape)
    R      = 5*np.ones(ids.shape)
    phi    = 45*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    theta  = np.random.rand(1,Nmrk)*2*np.pi
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    
    energy = Eprt*np.ones(ids.shape)
    pitch  = 0*np.ones(ids.shape)
    markers.write_hdf5_guidingcenters(fn[0], Nmrk, ids, mass, charge, 
                                      R, phi, z, energy, pitch, theta, weight, time)

    #subprocess.call(["cp", fn[0], fn[1]])
    #subprocess.call(["cp", fn[0], fn[2]])

    for i in range(0,3):
        # Options
        o = options.read_hdf5(fn[i])
        o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode[i]
        o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive[i]
        o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
        o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep[i]
        o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
        o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
        o["ENDCOND_SIMTIMELIM"]        = 0*o["ENDCOND_SIMTIMELIM"] + 1
        o["ENDCOND_MAX_SIM_TIME"]      = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime
        o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
        o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 0
        o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
        o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
        o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + simtime/1e5
        options.write_hdf5(fn[i],o)




def run():
    # Simulate. (The ascot5_momcoll is the main program "ascot5_main"
    # but compiled by setting A5_CCOL_NOENERGY to 1 in ascot5.h. This turns off
    # the energy part of the collisions.)
    subprocess.call(["./ascot5_main", "--in="+fn[0][0:-3]])
    #subprocess.call(["./ascot5_main", "--in="+fn[1][0:-3]])
    #subprocess.call(["./ascot5_main", "--in="+fn[2][0:-3]])


def check(plot=False):

    t, x, z = ([] for i in range(3))
    for i in range(0,len(fn)):
        if simmode[i] == 1:
            orb = ascot5.read_hdf5(fn[i],"orbits")["orbits"]["fo"]
        else:
            orb = ascot5.read_hdf5(fn[i],"orbits")["orbits"]["gc"]

        t.append(orb["time"])
        x.append( orb["R"] * np.cos(np.deg2rad(orb["phi"])) - 5 )
        #x.append(orb["R"])
        z.append(orb["z"])

    gamma = 1+Eprt*ELEMENTARY_CHARGE/(mprt*AMU2KG*constants.c*constants.c)
    vperp = np.sqrt(1-0.0*0.0) * np.sqrt(1-1/(gamma*gamma)) * constants.c
    
    #gyrolen = 0.9*np.sqrt(2*mprt*AMU2KG*Eprt*ELEMENTARY_CHARGE)/(Bnorm*ELEMENTARY_CHARGE*qprt)
    gyrolen = mprt*AMU2KG * gamma * vperp/(Bnorm*ELEMENTARY_CHARGE*np.abs(qprt))
    aphi = np.linspace(0,2*np.pi,1000)
    ax = gyrolen*np.cos(aphi)
    ay = gyrolen*np.sin(aphi)

    # Compare.

    # Plot if needed.
    if(plot):
        plt.figure()
        for i in range(0,1):
            #plt.plot(x[i],z[i],linestyle='.', linewidth=3,marker='.')
            plt.plot(x[i],z[i],linestyle='.', linewidth=3,marker='.')
            
        #plt.plot(ax,ay)
        plt.axis('equal')
        plt.show()


   

def clean():
    # Clean.
    for i in range(0,3):
        subprocess.call(["rm", "-f", fn[i]])

if __name__ == '__main__':
    init()
    run()
    check(plot=True)
    #clean()

