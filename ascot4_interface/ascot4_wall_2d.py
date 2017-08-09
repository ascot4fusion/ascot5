from pylab import *

import numpy as np

import sys
sys.path.append('../ui')
import ui_wall_2D as w2D

def read_wall_2d(fn):
    f = open(fn,'r')
    data = loadtxt(f,skiprows=1)
    f.close()
    str = {'r' : data[:,0], 'z' : data[:,1], 'divflag' : data[:,2]}
    return str

def write_wall_2d(f, w):
    w2D.write_hdf5(f, w['r'], w['z'])
