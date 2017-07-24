import numpy as np
import sys
sys.path.append('ui')
import ui_markers

def populate(fn, Rmin, Rmax, z0, N):
    ids    = np.arange(0,N) + 1
    r      = np.linspace(Rmin,Rmax,num=N)
    z      = np.ones(N)*z0
    phi    = np.ones(N)*0
    pitch  = np.ones(N)
    weight = np.ones(N)
    time   = np.ones(N)*0
    ui_markers.write_hdf5_fieldlines(fn, ids, r, phi, z, pitch, weight, time); 

if __name__ == "__main__":
    fn   = "ascot.h5"
    Rmin = 6.2
    Rmax = 8
    z0   = 0.6
    N    = 100
    populate(fn,Rmin,Rmax,z0,N)
