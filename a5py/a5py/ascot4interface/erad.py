from pylab import *
import numpy as np

def read_erad(fn):

    f = open(fn,'r')
    str = {'comm1' : f.readline()} #Skip comment line
    str['n_rho'] = f.readline()

    data = loadtxt(f)

    str['rho'] = data[:,0]
    str['dV_drho'] = data[:,1]
    # For data in format dV/rho, we can ignore effective minor radius
    
    return str
