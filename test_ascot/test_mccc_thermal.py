import h5py
import numpy as np
import subprocess
import sys 
from os import chdir
from os.path import dirname, realpath, sep, pardir
sys.path.append(dirname(realpath(__file__)) + sep + pardir + sep + 'ui')
import ui_B_TC
import ui_B_GS
import ui_E_TC
import ui_options
import ui_plasma_1D
import ui_wall_2D
import ui_markers
import ui_cleanout
import matplotlib.pyplot  as plt
import analyticBKGpsifun as psifun
from dists_coordinateTransform import vpavpe2Epitch
from ascot5_read import ascot5_read
import testfunctions as tf

def run():

    ############
    #INITIALIZE#
    ############

    fn = "ascot.h5"
    
    #Remove data from old runs
    ui_cleanout.clean(fn)

    #Create new options instance and set flags to zero:
    flags = ('ENABLE','ENDCOND')
    options = ui_options.Ui_optionsIO()
    ui_options.flagsToZero(options,flags)

    #Modify flags further as needed
    options.SIM_MODE = 1 #1 or 2, doesn't matter
    options.ENDCOND_ENERGYLIM = 1
    options.MIN_ENERGY = 10e2
    options.ENABLE_COULOMB_COLLISIONS = 1
    options.ENABLE_RZVparaVperp_DIST = 1

    options.DIST_RZVparaVperp_MIN_R = 1
    options.DIST_RZVparaVperp_MAX_R = 10
    options.DIST_RZVparaVperp_BIN_R = 1

    options.DIST_RZVparaVperp_MIN_Z = -10
    options.DIST_RZVparaVperp_MAX_Z = 10
    options.DIST_RZVparaVperp_BIN_Z = 1

    options.DIST_RZVparaVperp_MIN_VPARA = -1.5e7
    options.DIST_RZVparaVperp_MAX_VPARA = 1.5e7
    options.DIST_RZVparaVperp_BIN_VPARA = 40

    options.DIST_RZVparaVperp_MIN_VPERP = 0
    options.DIST_RZVparaVperp_MAX_VPERP = 1.5e7
    options.DIST_RZVparaVperp_BIN_VPERP = 20

    options.FIXEDSTEP_USERDEFINED = 1e-8
    options.FIXEDSTEP_USE_USERDEFINED = 1

    options.ENABLE_ORBIT_FOLLOWING = 0
    options.ENABLE_ORBITWRITE = 0


    axisr = 6.2;
    axisz = 0;
    psival = 0.5;
    rhoval = 0.5;
    B0 = np.array([1, 0, 0])
    B_dB = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
    E = np.array([0, 0, 0])
    Znum = np.array([1])
    Anum = np.array([1])
    rho = np.array([0, 0.5, 1, 2, 3])
    ndens = np.array([0, 0, 0, 0, 0])
    ntemp = np.array([0, 0, 0, 0, 0])
    edens = np.array([1, 1, 1, 1, 1]) * 1e20
    etemp = np.array([1, 1, 1, 1, 1]) * 1e0
    idens = np.array([1, 1, 1, 1, 1]) * 1e20
    itemp = np.array([1, 1, 1, 1, 1]) * 1e0
    Nrho = 5
    Nion = 1

    ndens = ndens.reshape(5,1)
    idens = idens.reshape(5,1)

    wr = np.array([0, 0, 8.3, 8.3, 0])
    wz = np.array([-4, 4, 4, -4, -4])

    ids    = np.array([1.0])
    mass   = np.array([4.002602])
    charge = np.array([2])
    weight = np.array([1])
    time   = np.array([0])

    #Particle
    rprt   = np.array([8.0123552])
    phiprt = np.array([352.52435])
    zprt   = np.array([0.0814756])
    vr     = np.array([8789194.5])
    vphi   = np.array([-3242767.0])
    vz     = np.array([-9101402.5])

    #GC
    r   = np.array([8.0123552])
    phi = np.array([352.52435])
    z   = np.array([0.0814756])
    pitch  = np.array([-0.7])
    energy = np.array([3.5e3])
    theta  = np.array([0.1])
    mlpitch = np.array([1.0])

    #B_GS
    axisr = 6.2
    axisz = 0
    B_phi0 = 5.3
    psimult = 200
    psicoef = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, -0.155])
    psi0 = psimult*psifun.psi0(axisr/axisr,axisz/axisr,psicoef[0],psicoef[1],psicoef[2],psicoef[3],psicoef[4],psicoef[5],psicoef[6],psicoef[7],psicoef[8],psicoef[9],psicoef[10],psicoef[11],psicoef[12])
    psi1 = 0


    #Input the new parameters in the hdf5 file
    ui_options.writeHdf5(options,fn)
    #ui_B_TC.write_hdf5(fn, B0, B_dB, axisr, axisz, psival, rhoval)
    #ui_B_GS.write_hdf5(fn, axisr, axisz, B_phi0, psi0, psi1, psimult, psicoef) 
    ui_B_GS.write_hdf5(fn, axisr, axisz, B_phi0, psi0, psi1, psimult, psicoef,Nripple=18,a0=2,delta0=0.5,alpha0=1.8) 
    #ui_B_GS.write_hdf5_B_3D(fn, axisr, axisz, B_phi0, psimult, psicoef,18, 2 ,3.8,0.01, np.array([3.9, 8.9, 100]), np.array([-5.0, 5.0, 200]), 360) 
    ui_E_TC.write_hdf5(fn, E) 
    ui_plasma_1D.write_hdf5(fn, Nrho, Nion, Znum, Anum, rho, ndens, ntemp, edens, etemp, idens, itemp)
    ui_wall_2D.write_hdf5(fn, wr, wz)
    #ui_markers.write_hdf5_particles(fn, ids, mass, charge, rprt, phiprt, zprt, vr, vphi, vz, weight, time);
    ui_markers.write_hdf5_guidingcenters(fn, ids+1, mass, charge, r, phi, z, energy, pitch, theta, weight, time);
    #ui_markers.write_hdf5_fieldlines(fn, ids+1, r, phi, z, mlpitch, weight, time);   

    #####
    #RUN#
    #####

    subprocess.call(['mv','ascot.h5','..'])
    chdir('..')
    subprocess.call(['make','ascot5_main','CC=mpicc','VERBOSE=1','NSIMD=1'])
    subprocess.call(['./ascot5_main'])

    #########
    #ANALYSE#
    #########

    f = h5py.File(fn, 'r')

    dist = ascot5_read(fn)['dists']
    mass = f['inistate/mass'][:]
    f.close()
    E_edges = np.linspace(0, 1e5, 41) 
    xi_edges = np.linspace(-1, 1, 21)
    E = np.linspace(0, 1e5, 40)

    Exidist = vpavpe2Epitch(dist,E_edges,xi_edges,mass)
    Exi = np.squeeze(Exidist['ordinate'])
    Edist = Exi.sum(axis=0)

    #print dist
    #print dist['ordinate']
    print Exidist['ordinate']
    #print Exidist['E']
    #print Exidist['xi']

    f = tf.MaxwellDistr(E)
    
    #plt.plot(E,f)
    plt.plot(Edist)
    plt.show()

if __name__ == '__main__':
    run()
