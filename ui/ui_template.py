import h5py
import numpy as np
import ui_B_TC
import ui_E_TC
import ui_options
import ui_plasma_1D
import ui_wall_2D
import ui_markers

fn = "ascot.h5"
f = h5py.File(fn, "a")
if "/options" in f:
    del f["options"]

if "/bfield" in f:
    del f["bfield"]

if "/efield" in f:
    del f["efield"]

if "/plasma" in f:
    del f["plasma"]

f.close();

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
anum   = np.array([4.0])
znum   = np.array([2.0])
weight = np.array([1.0])
rprt   = np.array([7.0123552])
phiprt = np.array([352.52435])
zprt   = np.array([0.0814756])
vr     = np.array([8789194.5])
vphi   = np.array([-3242767.0])
vz     = np.array([-9101402.5])

ui_options.writeHdf5(ui_options.ui_optionsIO(),fn)
ui_B_TC.write_hdf5(fn, B0, B_dB, axisr, axisz, psival, rhoval) 
ui_E_TC.write_hdf5(fn, E) 
ui_plasma_1D.write_hdf5(fn, Znum, Anum, rho, ndens, ntemp, edens, etemp, idens, itemp)
ui_wall_2D.write_hdf5(fn, wr, wz)
ui_markers.write_hdf5_particles(fn, ids, anum, znum, rprt, phiprt, zprt, vr, vphi, vz, weight);
