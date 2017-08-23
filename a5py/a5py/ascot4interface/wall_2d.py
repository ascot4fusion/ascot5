from pylab import *

import numpy as np

def read_wall_2d(fn):
    f = open(fn,'r')
    data = loadtxt(f,skiprows=1)
    f.close()
    str = {'r' : data[:,0], 'z' : data[:,1], 'divflag' : data[:,2]}
    return str
