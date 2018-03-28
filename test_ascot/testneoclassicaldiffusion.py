import numpy as np
import subprocess
import scipy.constants as constants
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.options as options
import a5py.ascot5io.B_GS as B_GS
import a5py.ascot5io.plasma_1D as plasma_1D
import a5py.ascot5io.markers as markers
import a5py.postprocessing.disttrasnformations as disttrans
import a5py.preprocessing.analyticequilibrium as psifun
import a5py.ascot5io.tools as tools
from time import sleep
from testcase import createbase

ELEMENTARY_CHARGE = constants.elementary_charge
AMU2KG = 1.66053904e-27

fn = ["neoclassGO.h5", "neoclassGCF.h5", "neoclassGCA.h5"]

scale = np.power(10,np.linspace(21.5,23.5,18))/1e21

simtime = 5e-3

mprt = 5.4858e-04

qprt = -1
Eprt = 1e3

Temp = 1e3
dens = 1e21

R0 = 6.2
z0 = 0
B_phi0 = 5.3
psi_mult = 200

# ITER-like equlibrium
#psi_coeff = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 
#                      3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 
#                      7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, 
#                      -0.155])

# ITER-like but circular equilibrium
psi_coeff = np.array([2.21808016e-02,  -1.28841781e-01,  -4.17718173e-02,
                      -6.22680280e-02,   6.20083978e-03,  -1.20524711e-03,
                      -3.70147050e-05,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      -0.155])
def init():
    
    # Test options
    simmode = [1, 2, 2]
    adaptive = [0, 0, 1]
    timestep = [1.e-9, 5e-8, 1e-8]
    tolorb = 1e-8
    tolcol = 1e-2
    Nmrk = 1000

    # Init testbase.
    Bxyz = np.array([0, 1, 0])
    Exyz = np.array([0, 0, 0])
    createbase(fn[0], Bxyz, Exyz, dens, Temp)

    

    # Init magnetic field.
    
    # Analytic equilibrium
    B_GS.write_hdf5(fn[0], R0, z0, B_phi0, psi_mult, psi_coeff)

    # Analytic equilibrium represented with 2D splines
    #B_GS.write_hdf5_B_2D(fn[0], R0, z0, B_phi0, psi_mult, psi_coeff, 
    #                     4, 8.5, 120, -4, 4, 200)

    # Analytic equilibrium represented with 3D splines
    #B_GS.write_hdf5_B_3D(fn[0], R0, z0, B_phi0, psi_mult, psi_coeff, 
    #                     18, 1, 1, delta0, 
    #                     4, 8.5, 120, -4, 4, 200, 0, 360, 360)

    # Markers
    ids    = np.linspace(1,Nmrk,Nmrk)
    mass   = mprt*np.ones(ids.shape)
    charge = qprt*np.ones(ids.shape)
    R      = 8*np.ones(ids.shape)
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
        for j in range(0,len(scale)):
            sleep(2)
            # Options
            qid,date = ascot5.get_qids(fn[i],"options")
            o = options.read_hdf5(fn[i],qid[0])
            o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode[i]
            o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive[i]
            o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
            o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep[i]/scale[j]
            o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
            o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
            o["ENDCOND_SIMTIMELIM"]        = 0*o["ENDCOND_SIMTIMELIM"] + 1
            o["ENDCOND_MAX_SIM_TIME"]      = 0*o["ENDCOND_MAX_SIM_TIME"] + simtime/scale[j]
            o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
            o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 1
            o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 0
            o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 0
            o["ORBITWRITE_INTERVAL"]       = 0*o["ORBITWRITE_INTERVAL"] + 1e-5
            options.write_hdf5(fn[i],o)

            # We wish to vary density, for which we need several plasma inputs
            Nrho = 3
            Nion = 1
            znum = np.array([1])
            anum = np.array([1])
            rho = np.array([0, 0.5, 1])
            ndens = 0*np.ones(rho.shape)
            ntemp = 0*np.ones(rho.shape)
            edens = np.ones(rho.shape)
            etemp = Temp*np.ones(rho.shape)
            idens = dens*np.ones((rho.size,1))
            itemp = Temp*np.ones(rho.shape)
            plasma_1D.write_hdf5(fn[i], Nrho, Nion, znum, anum, rho, ndens, ntemp, edens, etemp, idens*scale[j], itemp)


def run():
    # Simulate. (The ascot5_momcoll is the main program "ascot5_main"
    # but compiled by setting A5_CCOL_NOENERGY to 1 in ascot5.h. This turns off
    # the energy part of the collisions.)
    for i in range(1,2):
        qidpls,date = ascot5.get_qids(fn[i],"plasma")
        qidopt,date = ascot5.get_qids(fn[i],"options")
        for j in range(0,len(qidpls)-1):
            tools.setactivegroup(fn[i], "plasma", qidpls[j])
            tools.setactivegroup(fn[i], "options", qidopt[j])
            subprocess.call(["./ascot5_momcoll", "--in="+fn[i][0:-3]])
            sleep(2)


def check(plot=False):

    Rp = np.linspace(6.2,8.2,1000)
    zp = 0*np.ones(Rp.shape)
    c = psi_coeff
    psi = psi_mult*psifun.psi0(Rp/R0,zp/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])
    psi0 = psi_mult*psifun.psi0(1.0,0.0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])
    rhop = np.sqrt(np.absolute((psi-psi0)/psi0))

    DFO = []
    DGCF = []
    DGCA = []
    n = []

    DCL = []
    DB = []
    DPS = []
    
    for j in range(0,len(scale)):
        dr = []
        D  = []
        grid = np.linspace(-1,1,100)*5e-3
        for i in range(1,2):
            qid,date = ascot5.get_qids(fn[i],"results")
            state = ascot5.read_hdf5(fn[i],"states")["run-"+qid[j]]
            rf = np.interp(state["endstate"]["rho"],rhop,Rp)
            
            #state = ascot5.read_hdf5(fn[i],"states")["run-"+qid[j]]["inistate"]
            ri = np.interp(state["inistate"]["rho"],rhop,Rp)
            
            dr.append(np.histogram(ri-rf, bins=grid, density=True))
            D.append(np.var(ri-rf)/(simtime/scale[j]))
            
        grid = (grid[:-1] + grid[1:]) / 2
        
        clog = 14.7
        collfreq = 4.8e-14*dens*scale[j]*np.power(Temp,-3.0/2)*clog
        gyrolen = np.sqrt(2*mprt*AMU2KG*Eprt*ELEMENTARY_CHARGE)/(B_phi0*ELEMENTARY_CHARGE)
        
        eps = np.power(1.8/6.2,3.0/2.0)
        q = 2.2
        v = np.sqrt(2*Eprt*ELEMENTARY_CHARGE/(mprt*AMU2KG))
        
        #anD = (q/eps)*collfreq * np.power(gyrolen,2)
        #anD = (q*q/eps)*collfreq * np.power(gyrolen,2)
        #anD = (1+2*q*q)*collfreq * np.power(gyrolen,2)
    
        anD = q*np.power(gyrolen,2)*v / 6.2
        #anD = collfreq*6.2*q/v
        print((np.log10(dens*scale[j]),np.log10(anD*8*2),eps))
        
        anD = [0,0,0]
        anD[0] = collfreq * np.power(gyrolen,2)
        
        anD[1] = q*q*anD[0]*2*2*2*2
        anD[2] = np.power(6.2/1.9,3.0/2.0)*anD[1]
        
        # Compare.
        #print(anD)
        #print(D[0])
        #print(D[1])
        #print(D[2])

        n.append(dens*scale[j])
        #DFO.append(D[0])
        DGCF.append(D[0])
        #DGCA.append(D[2])
        DCL.append(anD[0])
        DPS.append(anD[1])
        DB.append(anD[2])
        
        # Plot if needed.
    if(plot):
        plt.figure()
        plt.plot(np.log10(n), np.log10(DGCF))
        #plt.plot(np.log10(n), np.log10(DCL))
        plt.plot(np.log10(n), np.log10(DB))
        plt.plot(np.log10(n), np.log10(DPS))
            
        plt.show()


    

def clean():
    # Clean.
    for i in range(0,3):
        subprocess.call(["rm", "-f", fn[i]])

if __name__ == '__main__':
    #init()
    #run()
    check(plot=True)
    #clean()

