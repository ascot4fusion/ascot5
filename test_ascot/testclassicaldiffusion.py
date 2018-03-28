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

fn = ["classGO.h5", "classGCF.h5", "classGCA.h5"]

ELEMENTARY_CHARGE = constants.elementary_charge
AMU2KG = 1.66053904e-27

simtime = 5e-6

mprt = 1.00727647
qprt = 1

Temp = 1e3
dens = 1e22
Bnorm = 5
Eprt = 1e3

def init():

    

    # Test options
    simmode = [1, 2, 2]
    adaptive = [0, 0, 1]
    timestep = [1.e-10, 1e-8, 1e-8]
    tolorb = 1e-8
    tolcol = 1e-1
    Nmrk = 1000



    # Init testbase.
    Bxyz = np.array([0, Bnorm, 0])
    Exyz = np.array([0, 0, 0])
    createbase(fn[0], Bxyz, Exyz, dens, Temp)

    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = mprt*np.ones(ids.shape)
    charge = qprt*np.ones(ids.shape)
    R      = 5*np.ones(ids.shape)
    phi    = 0*np.ones(ids.shape)
    z      = 0*np.ones(ids.shape)
    theta  = np.random.rand(1,Nmrk)*2*np.pi
    weight = 1*np.ones(ids.shape)
    time   = 0*np.ones(ids.shape)
    
    energy = Eprt*np.ones(ids.shape)
    pitch  = 1-2*np.random.rand(1,Nmrk)
    markers.write_hdf5_guidingcenters(fn[0], Nmrk, ids, mass, charge, 
                                      R, phi, z, energy, pitch, theta, weight, time)

    subprocess.call(["cp", fn[0], fn[1]])
    subprocess.call(["cp", fn[0], fn[2]])

    for i in range(0,3):
        # Options
        qid,date = ascot5.get_qids(fn[i],"options")
        o = options.read_hdf5(fn[i],qid[0])
        o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode[i]
        o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive[i]
        o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
        o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep[i]
        o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
        o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
        o["ENDCOND_SIMTIMELIM"]        = 0*o["ENDCOND_SIMTIMELIM"] + 1
        o["ENDCOND_MAX_SIM_TIME"]      = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime
        o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
        o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 1
        o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 1
        o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
        o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + 1e-5
        options.write_hdf5(fn[i],o)




def run():
    # Simulate. (The ascot5_momcoll is the main program "ascot5_main"
    # but compiled by setting A5_CCOL_NOENERGY to 1 in ascot5.h. This turns off
    # the energy part of the collisions.)
    subprocess.call(["./ascot5_main", "--in="+fn[0][0:-3]])
    subprocess.call(["./ascot5_momcoll", "--in="+fn[1][0:-3]])
    subprocess.call(["./ascot5_momcoll", "--in="+fn[2][0:-3]])


def check(plot=False):
    dr = []
    D  = []
    v = []

    grid = np.linspace(-1,1,50)*5e-3
    for i in range(0,3):
        qid,date = ascot5.get_qids(fn[i],"results")
        state = ascot5.read_hdf5(fn[i],"states")["run-"+qid[0]]["endstate"]
        xf = state["R"] * np.cos(np.deg2rad(state["phi"]))
        yf = state["R"] * np.sin(np.deg2rad(state["phi"]))
        zf = state["z"]
        rf = np.sqrt(np.power(xf,2) + np.power(zf,2))

        state = ascot5.read_hdf5(fn[i],"states")["run-"+qid[0]]["inistate"]
        xi = state["R"] * np.cos(np.deg2rad(state["phi"]))
        yi = state["R"] * np.sin(np.deg2rad(state["phi"]))
        zi = state["z"]
        ri = np.sqrt(np.power(xi,2) + np.power(zi,2))

        dr.append(np.histogram(zi-zf, bins=grid, density=True))
        D.append(np.var(ri-rf)/simtime)


    grid = (grid[:-1] + grid[1:]) / 2

    clog = 14.7
    collfreq = 4.8e-14*dens*np.power(Temp,-3.0/2)*clog
    gyrolen = np.sqrt(2*mprt*AMU2KG*Eprt*ELEMENTARY_CHARGE)/(Bnorm*ELEMENTARY_CHARGE)
    anD = collfreq * np.power(gyrolen,2)

    # Compare.
    print(anD)
    print(D[0])
    print(D[1])
    print(D[2])

    # Plot if needed.
    if(plot):
        plt.figure()
        for i in range(0,3):
            plt.plot(grid,dr[i][0])
            
        plt.show()


   

def clean():
    # Clean.
    for i in range(0,3):
        subprocess.call(["rm", "-f", fn[i]])

if __name__ == '__main__':
    init()
    run()
    check(plot=True)
    clean()

