import numpy as np

def read_wall_2d(fn):
    with open(fn,'r') as f:
        h5data = np.loadtxt(f,skiprows=1)
    data = {'r' : h5data[:,0], 'z' : h5data[:,1], 'divflag' : h5data[:,2]}
    return data
