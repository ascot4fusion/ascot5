import numpy as np
import sys
sys.path.append('ui')
import ui_markers

def populate(fn, Rmin, Rmax, z0, N):
    N = (92-10*2)*(168-10*2)
    ids    = np.arange(0,N) + 1
    r      = np.linspace(3.52,8.98,num=(92))
    z      = np.linspace(-5.01,5.01,num=(168))
    xv, yv = np.meshgrid(r[10:-10]-0.011, z[10:-10]-0.011)
    phi    = np.ones(N)*7
    pitch  = np.ones(N)
    weight = np.ones(N)
    time   = np.ones(N)*0
    ui_markers.write_hdf5_fieldlines(fn, ids, xv.flatten(), phi, yv.flatten(), pitch, weight, time); 

if __name__ == "__main__":
    fn   = "ascot.h5"
    Rmin = 6.2
    Rmax = 8
    z0   = 0.6
    N    = 100
    populate(fn,Rmin,Rmax,z0,N)
