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
import testfunctions as tf
import scipy.constants as const

#Invariants
def run(B2D):
    ############
    #INITIALIZE#
    ############

    fn = 'ascot.h5'

    #Remove data from old runs
    ui_cleanout.clean(fn)

    #Create new options instance and set flags to zero:
    flags = ('ENABLE','ENDCOND')
    options = ui_options.Ui_optionsIO()
    ui_options.flagsToZero(options,flags)

    #Modify flags further as needed
    options.SIM_MODE               = 2
    options.ENDCOND_SIMTIMELIM     = 1
    options.ENDCOND_MAX_SIM_TIME   = 1e-3
    options.ENABLE_ORBIT_FOLLOWING = 1
    options.ENABLE_ORBITWRITE      = 1
    options.ORBITWRITE_INTERVAL    = 1e-1
    options.ORBITWRITE_MODE           = 1
    options.FIXEDSTEP_USE_USERDEFINED = 1
    options.FIXEDSTEP_USERDEFINED = 1e-8

    axisr = 0
    axisz = 0
    psival = 0.5
    rhoval = 0.5
    B0 = np.array([1, 0, 0])
    B_dB = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
    E = np.array([0, 0, 0])
    Znum = np.array([1])
    Anum = np.array([1])
    rho = np.array([0, 0.5, 1, 2, 3])
    ndens = np.array([0, 0, 0, 0, 0])
    ntemp = np.array([0, 0, 0, 0, 0])
    edens = np.array([1, 1, 1, 1, 1]) * 1e20
    etemp = np.array([1, 1, 1, 1, 1]) * 1e4
    idens = np.array([1, 1, 1, 1, 1]) * 1e20
    itemp = np.array([1, 1, 1, 1, 1]) * 1e4
    Nrho = 5
    Nion = 1

    ndens = ndens.reshape(5,1)
    idens = idens.reshape(5,1)

    wr = np.array([0, 0, 100, 100, 0])
    wz = np.array([-50, 50, 50, -50, -50])

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
    pitch  = np.array([0.7])
    energy = np.array([3.5e6])
    theta  = np.array([0.1])
    mlpitch = np.array([1.0])

    #B_GS
    axisr = 6.2
    axisz = 0
    psi0 = -0.0365
    psi1 = 0
    B_phi0 = 5.3
    psimult = 200
    psicoef = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, -0.155])

    #Write variables in the hdf5 file
    ui_options.writeHdf5(options,fn)
    #ui_B_TC.write_hdf5(fn, B0, B_dB, axisr, axisz, psival, rhoval) 
    if B2D:
        ui_B_GS.write_hdf5_B_2D(fn, axisr, axisz, B_phi0, psimult, psicoef, np.array([3.9, 8.9, 400]), np.array([-5.0, 5.0, 800])) 
    else:
        ui_B_GS.write_hdf5(fn, axisr, axisz, B_phi0, psi0, psi1, psimult, psicoef) 
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

    #Inistate
    mass_ini = tf.amu2kg(f['inistate/mass'][:])
    charge_ini = tf.e2C(f['inistate/charge'][:])

    R_ini = f['inistate/Rprt'][:]
    phi_ini = f['inistate/phi'][:]
    z_ini = f['inistate/zprt'][:]

    vR_ini = f['inistate/vR'][:]
    vphi_ini = f['inistate/vphi'][:] 
    vxy_ini = tf.v(vR_ini,vphi_ini,phi_ini)
    vz_ini = f['inistate/vz'][:]
    vv_ini = [vxy_ini[0][0],vxy_ini[1][0],vz_ini[0]]
    v_tot_ini = tf.v3D(vxy_ini,vz_ini)
    
    BR_ini = f['inistate/B_R'][:]
    Bphi_ini = f['inistate/B_phi'][:]
    Bxy_ini = tf.v(BR_ini,Bphi_ini,phi_ini)
    Bz_ini = f['inistate/B_z'][0]
    B_ini = [Bxy_ini[0][0],Bxy_ini[1][0],Bz_ini]
    B_tot_ini = tf.v3D(Bxy_ini,Bz_ini)
    
    v_ini = tf.vectProjection(vv_ini,B_ini) # = (v_para, v_perp)
    v_para_ini = tf.v3D((v_ini[0][0],v_ini[0][1]),v_ini[0][2])
    v_perp_ini  =tf.v3D((v_ini[1][0],v_ini[1][1]),v_ini[1][2])
    
    #Endstate
    mass_end = tf.amu2kg(f['endstate/mass'][:])
    charge_end = tf.e2C(f['endstate/charge'][:])

    R_end = f['endstate/Rprt'][:]
    phi_end = f['endstate/phiprt'][:]
    z_end = f['endstate/zprt'][:]

    vR_end = f['endstate/vR'][:]
    vphi_end = f['endstate/vphi'][:] 
    vxy_end = tf.v(vR_end,vphi_end,phi_end)
    vz_end = f['endstate/vz'][:]
    vv_end = [vxy_end[0][0],vxy_end[1][0],vz_end[0]]
    v_tot_end = tf.v3D(vxy_end,vz_end)
    
    BR_end = f['endstate/B_R'][:]
    Bphi_end = f['endstate/B_phi'][:]
    Bxy_end = tf.v(BR_end,Bphi_end,phi_end)
    Bz_end = f['endstate/B_z'][0]
    B_end = [Bxy_end[0][0],Bxy_end[1][0],Bz_end]
    B_tot_end = tf.v3D(Bxy_end,Bz_end)
    
    v_end = tf.vectProjection(vv_end,B_end) # = (v_para, v_perp)
    v_para_end = tf.v3D((v_end[0][0],v_end[0][1]),v_end[0][2])
    v_perp_end  =tf.v3D((v_end[1][0],v_end[1][1]),v_end[1][2])
    
    #Energy
    energy_ini = tf.totalEnergy(mass_ini,v_para_ini,v_perp_ini)
    energy_end = tf.totalEnergy(mass_end,v_para_end,v_perp_end)

    #Magnetic moment
    mu_ini = f['inistate/mu'][:]
    mu_end = f['endstate/mu'][:]

    #Canonical momentum
    psi_mult = f['bfield/B_GS/psi_mult'][()]
    c = f['bfield/B_GS/psi_coeff'][:]
    R0 = f['bfield/B_GS/R0'][()]
    z0 = f['bfield/B_GS/z0'][()]
    psi = psi_mult*psifun.psi0(R0/R0,z0/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])
    pphi_ini = tf.canonicalMomentum(mass_ini,R_ini,v_para_ini,charge_ini,psi)
    pphi_end = tf.canonicalMomentum(mass_end,R_end,v_para_end,charge_end,psi)

    print 'Initial total energy: %s' % energy_ini
    print 'End total energy: %s' % energy_end
    print 'Initial magnetic moment: %s' % mu_ini
    print 'End magnetic moment: %s' % mu_end
    print 'Initial toroidal canonical momentum: %s' % pphi_ini
    print 'End toroidal canonical momentum: %s' % pphi_end
    
    f.close()

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'True':
        run(True)
    else:
        run(False)
