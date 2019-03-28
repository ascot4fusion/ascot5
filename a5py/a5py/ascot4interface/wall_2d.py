import numpy as np

def read_wall_2d(fn):
    with open(fn,'r') as f:
        data = np.loadtxt(f,skiprows=1)
    str = {'r' : data[:,0], 'z' : data[:,1], 'divflag' : data[:,2]}
    return str
