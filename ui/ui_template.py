import h5py
import numpy as np
import ui_B_TC
import ui_B_GS
import ui_E_TC
import ui_options
import ui_plasma_1D
import ui_wall_2D
import ui_markers
import ui_cleanout
import subprocess

fn = "ascot.h5"

#Remove data from old runs
ui_cleanout.clean(fn)

#Create new options instance and set flags to zero:
flags = ('ENABLE','ENDCOND')
options = ui_options.Ui_optionsIO()
ui_options.flagsToZero(options,flags)

#Modify flags further as needed

axisr = 0;
axisz = 0;
psival = 0.5;
rhoval = 0.5;
B0 = np.array([1, 0, 0])
B_dB = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
E = np.array([1, 0, 0])
Znum = 1
Anum = 1
rho = np.array([0, 0.5, 1, 2, 3])
ndens = np.array([0, 0, 0, 0, 0])
ntemp = np.array([0, 0, 0, 0, 0])
edens = np.array([1, 1, 1, 1, 1]) * 1e20
etemp = np.array([1, 1, 1, 1, 1]) * 1e4
idens = np.array([1, 1, 1, 1, 1]) * 1e20
itemp = np.array([1, 1, 1, 1, 1]) * 1e4

ndens = ndens.reshape(5,1)
idens = idens.reshape(5,1)

wr = np.array([0, 0, 100, 100, 0])
wz = np.array([-50, 50, 50, -50, -50])

ids    = np.array([1.0])
mass   = np.array([4.002602])
charge = np.array([2])
weight = np.array([1])
time   = np.array([0])

rprt   = np.array([7.0123552])
phiprt = np.array([352.52435])
zprt   = np.array([0.0814756])
vr     = np.array([8789194.5])
vphi   = np.array([-3242767.0])
vz     = np.array([-9101402.5])

r   = np.array([7.0123552])
phi = np.array([352.52435])
z   = np.array([0.0814756])
pitch  = np.array([0.7])
energy = np.array([3.5e6])
theta  = np.array([0.1])
mlpitch = np.array([1.0])

axisr = 6.2
axisz = 0
psi0 = -0.0365
psi1 = 0
B_phi0 = 5.3
psimult = 200
psicoef = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, -0.155])

#Input the new parameters in the hdf5 file
ui_options.writeHdf5(ui_options.Ui_optionsIO(),fn)
#ui_B_TC.write_hdf5(fn, B0, B_dB, axisr, axisz, psival, rhoval)
#ui_B_GS.write_hdf5(fn, axisr, axisz, B_phi0, psi0, psi1, psimult, psicoef) 
ui_B_GS.write_hdf5_B_2D(fn, axisr, axisz, B_phi0, psimult, psicoef, np.array([3.9, 8.9, 400]), np.array([-5.0, 5.0, 800])) 
ui_E_TC.write_hdf5(fn, E) 
ui_plasma_1D.write_hdf5(fn, Znum, Anum, rho, ndens, ntemp, edens, etemp, idens, itemp)
ui_wall_2D.write_hdf5(fn, wr, wz)
#ui_markers.write_hdf5_particles(fn, ids, mass, charge, rprt, phiprt, zprt, vr, vphi, vz, weight, time);
#ui_markers.write_hdf5_guidingcenters(fn, ids+1, mass, charge, r, phi, z, energy, pitch, theta, weight, time);
#ui_markers.write_hdf5_fieldlines(fn, ids+1, r, phi, z, mlpitch, weight, time);   

#Compile and run ASCOT5
subprocess.call(['make','ascot5_main','CC=mpicc','VERBOSE=1','NSIMD=1'])
subprocess.call(['./ascot5_main'])

#Analyse the data returned by ASCOT5 
f = h5py.File(fn, 'r')
f.close()
