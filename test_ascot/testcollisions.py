import numpy as np
import subprocess
import scipy.constants as constants
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.options as options
import a5py.ascot5io.B_GS as B_GS
import a5py.ascot5io.markers as markers
import a5py.postprocessing.disttrasnformations as disttrans
from a5py.postprocessing.state import gather
from testcase import createbase

ELEMENTARY_CHARGE = constants.elementary_charge
AMU2KG = 1.66053904e-27

dens = 1e20
Temp = 1e3
mprt = 1.00727647

E0 = 1e5

fn = ["coll_thGO.h5", "coll_thGCF.h5", "coll_thGCA.h5", "coll_sdGO.h5", "coll_sdGCF.h5", "coll_sdGCA.h5"]

def init(fast=False):

    # Test options
    simmode = [1, 2, 2, 1, 2, 2]  
    adaptive = [0, 0, 1, 0, 0, 1]
    timestep = [1e-6, 1e-6, 1e-6, 1e-7, 1e-7, 1e-7]
    tolorb = 1e-8
    tolcol = 2e-1
    orbfollowing = 0
    Nmrksd = 1000
    Nmrkth = 2000

    # Proton
    m = mprt
    q = 1

    # Init testbase.

    Bxyz = np.array([1, 0, 0])
    Exyz = np.array([0, 0, 0])
    createbase(fn[0], Bxyz, Exyz, dens, Temp)

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

    # Clone test case base
    for i in range(1,6):
        subprocess.call(["cp",fn[0],fn[i]])

    # Common options and markers
    for i in range(0,6):
        qid,date = ascot5.get_qids(fn[0],"options")
        o = options.read_hdf5(fn[0],qid[0])
        o["SIM_MODE"]                  = 0*o["SIM_MODE"] + simmode[i]
        o["ENABLE_ADAPTIVE"]           = 0*o["ENABLE_ADAPTIVE"] + adaptive[i]
        o["FIXEDSTEP_USE_USERDEFINED"] = 0*o["FIXEDSTEP_USE_USERDEFINED"] + 1
        o["FIXEDSTEP_USERDEFINED"]     = 0*o["FIXEDSTEP_USERDEFINED"] + timestep[i]
        o["ADAPTIVE_TOL_ORBIT"]        = 0*o["ADAPTIVE_TOL_ORBIT"] + tolorb
        o["ADAPTIVE_TOL_CCOL"]         = 0*o["ADAPTIVE_TOL_CCOL"] + tolcol
        o["ENABLE_ORBIT_FOLLOWING"]    = 0*o["ENABLE_ORBIT_FOLLOWING"] + orbfollowing
        o["ENABLE_COULOMB_COLLISIONS"] = 0*o["ENABLE_COULOMB_COLLISIONS"] + 1
        
        o["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"]  = 0*o["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"] + 1
        o["DIST_MIN_R"]   = 0*o["DIST_MIN_R"] + 4
        o["DIST_MAX_R"]   = 0*o["DIST_MAX_R"] + 8
        o["DIST_NBIN_R"]  = 0*o["DIST_NBIN_R"] + 1
        o["DIST_MIN_z"]   = 0*o["DIST_MIN_z"] - 5
        o["DIST_MAX_z"]   = 0*o["DIST_MAX_z"] + 5
        o["DIST_NBIN_z"]  = 0*o["DIST_NBIN_z"] + 1
        o["ENABLE_ORBITWRITE"]         = 0*o["ENABLE_ORBITWRITE"] + 0
        o["ORBITWRITE_MODE"]           = 0*o["ORBITWRITE_MODE"] + 1
        options.write_hdf5(fn[i],o)

    # Thermal options and markers
    for i in range(0,3):
        qid,date = ascot5.get_qids(fn[i],"options")
        o = options.read_hdf5(fn[i],qid[0])
        o["ENDCOND_SIMTIMELIM"]          = 0*o["ENDCOND_SIMTIMELIM"] + 1
        o["ENDCOND_MAX_SIM_TIME"]        = 0*o["ENDCOND_MAX_SIM_TIME"] + 1e-3
        o["ORBITWRITE_INTERVAL"]         = 0*o["ORBITWRITE_INTERVAL"] + 1e-5
        o["DIST_MIN_vpa"]  = 0*o["DIST_MIN_vpa"] - 2e6
        o["DIST_MAX_vpa"]  = 0*o["DIST_MAX_vpa"] + 2e6
        o["DIST_NBIN_vpa"] = 0*o["DIST_NBIN_vpa"] + 100
        o["DIST_MIN_vpe"]  = 0*o["DIST_MIN_vpe"] + 0
        o["DIST_MAX_vpe"]  = 0*o["DIST_MAX_vpe"] + 2e6
        o["DIST_NBIN_vpe"] = 0*o["DIST_NBIN_vpe"] + 50
        options.write_hdf5(fn[i],o)

        ids    = np.linspace(1,Nmrkth,Nmrkth)
        mass   = m*np.ones(ids.shape)
        charge = q*np.ones(ids.shape)
        R      = 6.4*np.ones(ids.shape)
        phi    = 0*np.ones(ids.shape)
        z      = 0*np.ones(ids.shape)
        theta  = 0*np.ones(ids.shape)
        weight = 1*np.ones(ids.shape)
        time   = 0*np.ones(ids.shape)

        energy = 1.0e3*np.ones(ids.shape)
        pitch  = 0.9*np.ones(ids.shape)
        markers.write_hdf5_guidingcenters(fn[i], Nmrkth, ids, mass, charge, 
                                          R, phi, z, energy, pitch, theta, weight, time)

    # Slowing down options and markers
    for i in range(3,6):
        qid,date = ascot5.get_qids(fn[i],"options")
        o = options.read_hdf5(fn[i],qid[0])
        o["ENDCOND_ENERGYLIM"]                = 0*o["ENDCOND_ENERGYLIM"] + 1
        o["ENDCOND_MIN_ENERGY"]               = 0*o["ENDCOND_MIN_ENERGY"] + 1e3
        o["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = 0*o["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] + 0.01
        o["ORBITWRITE_INTERVAL"]              = 0*o["ORBITWRITE_INTERVAL"] + 1e-4
        o["DIST_MIN_vpa"]      = 0*o["DIST_MIN_vpa"] - 5e6
        o["DIST_MAX_vpa"]      = 0*o["DIST_MAX_vpa"] + 5e6
        o["DIST_NBIN_vpa"]     = 0*o["DIST_NBIN_vpa"] + 200
        o["DIST_MIN_vpe"]      = 0*o["DIST_MIN_vpe"] + 0
        o["DIST_MAX_vpe"]      = 0*o["DIST_MAX_vpe"] + 5e6
        o["DIST_NBIN_vpe"]     = 0*o["DIST_NBIN_vpe"] + 100
    
        options.write_hdf5(fn[i],o)

        ids    = np.linspace(1,Nmrksd,Nmrksd)
        mass   = m*np.ones(ids.shape)
        charge = q*np.ones(ids.shape)
        R      = 6.4*np.ones(ids.shape)
        phi    = 0*np.ones(ids.shape)
        z      = 0*np.ones(ids.shape)
        theta  = 0*np.ones(ids.shape)
        weight = 1*np.ones(ids.shape)
        time   = 0*np.ones(ids.shape)
    
        energy = E0*np.ones(ids.shape)
        pitch  = 1-2*np.random.rand(Nmrksd)
        markers.write_hdf5_guidingcenters(fn[i], Nmrksd, ids, mass, charge, 
                                          R, phi, z, energy, pitch, theta, weight, time)



def run():
    # Simulate.
    for i in range(0,3):
        subprocess.call(["./ascot5_main", "--in="+fn[i][0:-3]])
    for i in range(3,6):
        subprocess.call(["./ascot5_main", "--in="+fn[i][0:-3]])

def check(plot=False):

    thEgrid    = np.linspace(0,6*1e3,50)
    sdEgrid    = np.linspace(0,1.1*1e5,50)
    sdtimegrid = np.linspace(0,0.03,50)
    pitchgrid  = np.linspace(-1,1,50)

    thEkin, thpitch, sdtime, sdpitch, sdvdist, sdEdist = ([] for i in range(6))

    # Gather thermal process data
    for i in range(0,3):
        qid,date = ascot5.get_qids(fn[i],"results")
        s = ascot5.read_hdf5(fn[i],"states")["run-"+qid[0]]["endstate"]
        gather(s)

        thEkin.append(np.histogram(s["Ekin"],   bins=thEgrid, density=True))
        thpitch.append(np.histogram(s["pitch"], bins=pitchgrid, density=True))

    for i in range(3,6):
        qid,date = ascot5.get_qids(fn[i],"results")
        s = ascot5.read_hdf5(fn[i],"states")["run-"+qid[0]]["endstate"]
        gather(s)
        
        sdpitch.append(np.histogram(s["pitch"], bins=pitchgrid, density=True))
        sdtime.append(np.histogram(s["time"], bins=sdtimegrid, density=True))
        
        dist = ascot5.read_hdf5(fn[i],"dists")["run-"+qid[0]]["dists"]["R_phi_z_vpa_vpe_t_q"]
        vpagrid = dist["vpa"]
        vpegrid = dist["vpe"]
        sdvdist.append(np.transpose(np.squeeze(dist["ordinate"])))

        RzExi = disttrans.vpavpe2Epitch(dist, sdEgrid, pitchgrid, mprt)
        RzExi = (np.squeeze(RzExi['ordinate'])).sum(axis=(1))
        sdEdist.append(RzExi/np.sum( RzExi*(sdEgrid[1]-sdEgrid[0]) ))

    thEgrid    = (thEgrid[:-1] + thEgrid[1:]) / 2
    pitchgrid  = (pitchgrid[:-1] + pitchgrid[1:]) / 2
    sdEgrid    = (sdEgrid[:-1] + sdEgrid[1:]) / 2
    sdtimegrid = (sdtimegrid[:-1] + sdtimegrid[1:]) / 2

    anthEkin = np.sqrt(thEgrid)*np.exp(-thEgrid/Temp)
    anthEkin = anthEkin/np.sum(anthEkin*(thEgrid[1]-thEgrid[0]))

    anthpitch = np.ones(pitchgrid.shape)
    anthpitch = anthpitch/np.sum(anthpitch*(pitchgrid[1]-pitchgrid[0]))
    ansdpitch = anthpitch

    v0 = E2v(E0*ELEMENTARY_CHARGE, mprt*AMU2KG)
    v = E2v(sdEgrid*ELEMENTARY_CHARGE, mprt*AMU2KG)
    vth = E2v(Temp*ELEMENTARY_CHARGE, constants.electron_mass)
    vcrit = vth * np.power( (3.0*np.sqrt(np.pi)/4.0) * constants.electron_mass / constants.proton_mass , 1/3.0)
    ansdEkin = heaviside(v0-v, 1) *  v / (np.power(vcrit,3) + np.power(v,3))
    ansdEkin = ansdEkin/np.sum(ansdEkin*(sdEgrid[1]-sdEgrid[0]))

    # Compare.


    # Plot if needed.
    if(plot):
        c = ['b', 'g', 'r']
        f  = plt.figure()
        a1 = f.add_subplot(3,3,1)
        a2 = f.add_subplot(3,3,2)
        a3 = f.add_subplot(3,3,3)
        a4 = f.add_subplot(3,3,4)
        a5 = f.add_subplot(3,3,5)
        a6 = f.add_subplot(3,3,6)
        a7 = f.add_subplot(3,3,7)
        a8 = f.add_subplot(3,3,8)
        a9 = f.add_subplot(3,3,9)

        a1.plot(thEgrid,anthEkin)
        a2.plot(pitchgrid,anthpitch)

        a3.plot(sdEgrid,ansdEkin)
        a4.plot(0,0)
        a5.plot(pitchgrid,ansdpitch)

        for i in range(0,3):
            a1.plot(thEgrid,thEkin[i][0])
            a2.plot(pitchgrid,thpitch[i][0])

        for i in range(0,3):
            a3.plot(sdEgrid,sdEdist[i])
            a4.plot(sdtimegrid,sdtime[i][0])
            a5.plot(pitchgrid,sdpitch[i][0])

            
        a1.set_xlabel('Energy (eV)'), a1.set_ylabel('f(E) (a.u.)'), a1.set_title('Thermal final energy')
        a2.set_xlabel('pitch (1)'), a2.set_ylabel('f(xi) (a.u.)'), a2.set_title('Thermal final pitch')
        a3.set_xlabel('Energy (eV)'), a3.set_ylabel('f(E) (a.u.)'), a3.set_title('Slowing-down energy dist')
        a5.set_xlabel('pitch (1)'), a5.set_ylabel('f(xi) (a.u.)'), a5.set_title('Slowing-down final pitch')

        a4.set_xlabel('Time (s)'), a4.set_ylabel('f(t) (a.u.)'), a4.set_title('Slowing-down time')
        

        a7.pcolor(vpagrid,vpegrid,sdvdist[0]), a7.set_xlabel('vpa (m/s)'), a7.set_ylabel('vpe (m/s)'), a7.set_title('GO: Slowing down dist')
        a8.pcolor(vpagrid,vpegrid,sdvdist[1]), a8.set_xlabel('vpa (m/s)'), a8.set_ylabel('vpe (m/s)'), a7.set_title('GCF: Slowing down dist')
        a9.pcolor(vpagrid,vpegrid,sdvdist[2]), a9.set_xlabel('vpa (m/s)'), a9.set_ylabel('vpe (m/s)'), a7.set_title('GCA: Slowing down dist')

        l1, = a6.plot(0, 0, label='Analytic')
        l2, = a6.plot(0, 0, label='GO')
        l3, = a6.plot(0, 0, label='GC Fixed')
        l4, = a6.plot(0, 0, label='GC Adaptive')
        handles= [l1, l2, l3, l4]
        labels = ['Analytic', 'GO', 'GC Fixed', 'GC Adaptive']
        a6.legend(handles, labels)
        a6.get_xaxis().set_ticks([])
        a6.get_yaxis().set_ticks([])

        plt.show()

    

def clean():
    # Clean.
    for i in range(0,6):
        subprocess.call(["rm", "-f", fn[i]])

def E2v(E,m):
    gamma = 1 + E/(m*constants.c*constants.c)
    v = np.sqrt(1 - 1.0/(gamma*gamma))*constants.c
    return v

def heaviside(x,x0):
    y = np.zeros(x.shape)
    y[x>=x0] = 1.0
    return y

if __name__ == '__main__':
    print("\n\
    This test validates the Coulomb collision operator.\n\
    Markers are simulated in ITER-like background with collisions,\n\
    and orbit-following if the slow-simulation mode is turned on.\n\
    \n\
    Two tests are performed and repeated with GO, GC fixed and GC adaptive modes.\n\
    Thermal equilibrium test checks that ions relax to Maxwellian distribution.\n\
    Slowing-down simulation checks that correct slowing-down distribution is produced.\n\
    ")

    init()
    run()
    check(True)
    clean()
