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
import testfunctions as tf

#B=1, E=5e5
def run():
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
    options.SIM_MODE               = 1
    options.ENDCOND_SIMTIMELIM     = 1
    options.ENDCOND_MAX_SIM_TIME   = 2e-7
    options.ENABLE_ORBIT_FOLLOWING = 1
    options.ENABLE_ORBITWRITE      = 1
    options.ORBITWRITE_INTERVAL    = 1.0e-11
    options.ORBITWRITE_MODE           = 1
    options.FIXEDSTEP_USE_USERDEFINED = 1
    options.FIXEDSTEP_USERDEFINED = 1e-11

    #Variables
    axisr = 0
    axisz = 0
    psival = 0.5
    rhoval = 0.5
    B0 = np.array([0, 0, 1.0])
    B_dB = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
    E = np.array([0, 5.0e5, 0])
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
    mass   = np.array([1.0])
    charge = np.array([1.0])
    weight = np.array([1.0])
    time = 0

    #Particle
    rprt   = np.array([7.0123552])
    phiprt = np.array([352.52435])
    zprt   = np.array([0.0814756])
    vr     = np.array([8789194.5])
    vphi   = np.array([-3242767.0])
    vz     = np.array([-9101402.5])

    #GC
    r       = np.array([7.0123552])
    phi     = np.array([352.52435])
    z       = np.array([0.0814756])
    pitch   = np.array([0.7])
    energy  = np.array([3.5e6])
    theta   = np.array([0.1])
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
    ui_B_TC.write_hdf5(fn, B0, B_dB, axisr, axisz, psival, rhoval) 
    #ui_B_GS.write_hdf5(fn, axisr, axisz, B_phi0, psi0, psi1, psimult, psicoef) 
    #ui_B_GS.write_hdf5_B_2D(fn, axisr, axisz, B_phi0, psimult, psicoef, np.array([3.9, 8.9, 400]), np.array([-5.0, 5.0, 800])) 
    ui_E_TC.write_hdf5(fn, E) 
    ui_plasma_1D.write_hdf5(fn, Nrho, Nion, Znum, Anum, rho, ndens, ntemp, edens, etemp, idens, itemp)
    ui_wall_2D.write_hdf5(fn, wr, wz)
    ui_markers.write_hdf5_particles(fn, ids, mass, charge, rprt, phiprt, zprt, vr, vphi, vz, weight, time)
    #ui_markers.write_hdf5_guidingcenters(fn, ids, mass, charge, r, phi, z, energy, pitch, theta, weight, time)
    #ui_markers.write_hdf5_fieldlines(fn, ids, r, phi, z, mlpitch, weight, time)

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

    massSI = tf.amu2kg(mass)
    chargeSI = tf.e2C(charge)

    mass_out = f['orbits/fo/mass'][:]
    charge_out = f['orbits/fo/charge'][:]
    time_out = f['orbits/fo/time'][:]

    #Particle coordinates in every timestep:
    r_out = f['orbits/fo/R'][:]
    phi_out = f['orbits/fo/phi'][:]
    phi_rad = np.deg2rad(phi_out)
    x_out = tf.xFromPol(r_out,phi_rad)
    y_out = tf.yFromPol(r_out,phi_rad)
    z_out = f['orbits/fo/z'][:]

    #Particle velocity in every timestep:
    vr_out = f['orbits/fo/v_R'][:]
    vphi_out = f['orbits/fo/v_phi'][:]
    vphi_rad = np.deg2rad(vphi_out)
    vz_out = f['orbits/fo/v_z'][:]

    v = tf.v(vr_out, vphi_rad, phi_rad)
    v_xy = tf.v2D(v)

    #initial GC coordinates:
    rGCini = f['inistate/R'][:]
    phiGCini = f['inistate/phi'][:]
    phiGCini_rad = np.deg2rad(phiGCini)
    xGCini = tf.xFromPol(rGCini,phiGCini_rad)
    yGCini = tf.yFromPol(rGCini,phiGCini_rad)
    zGCini = f['inistate/z'][:]

    #end GC coordinates:
    rGCend = f['endstate/R'][:]
    phiGCend = f['endstate/phi'][:]
    phiGCend_rad = np.deg2rad(phiGCend)
    xGCend = tf.xFromPol(rGCend,phiGCend_rad)
    yGCend = tf.yFromPol(rGCend,phiGCend_rad)
    zGCend = f['endstate/z'][:]

    f.close()

    #Analytic larmor radius
    radius_a =  tf.Larmor_a(massSI,v_xy,chargeSI,B0[2])[0]

    #Analytic period of rotation
    T_a = tf.T_a(radius_a,massSI*v_xy,massSI)[0]

    #Analytic drift velocity of GC
    #vGC_a = tf.vGC_a(vz_out[0],B0,E)
    #vGC_a = tf.vGC_EcrossB_a(E,B0)
    vCG_a = tf.vGC_a(vz_out,B0,E)

    #Actual drift velocity of GC
    vGC = tf.vGClarmor(x_out,y_out,vr_out,vphi_rad,phi_rad,charge,radius_a, time_out)
    #vGC = tf.vGC(xGCini,xGCend,time_out)

    #Actual period of rotation
    T = tf.T(v[0],time_out)

    #Actual larmor radius
    #radius = tf.Larmor2(y_out, v[0])
    radius = tf.Larmor(x_out, y_out, T_a, time_out)

    print 'Analytic larmor radius: %s' % radius_a
    print 'Actual larmor radius: %s' % radius
    print 'Analytic period of rotation: %s' % T_a
    print 'Actual period of rotation: %s' % T
    print 'Analytic drift velocity of GC: %s' % vCG_a
    print 'Actual drift velocity of GC %s' % vGC

if __name__ == '__main__':
    run()

