from pylab import *

def read_wall_2d(fn):
    f = open(fn,'r')
    data = loadtxt(f,skiprows=1)
    f.close()
    str = {'r' : data[:,0], 'z' : data[:,1], 'divflag' : data[:,2]}
    return str

def write_wall_2d(f, w):

    f.create_group('wall/2D')
    f['wall/2D'].attrs['n'] = w['r'].size
    f.create_dataset('wall/2D/r', data=w['r'])
    f.create_dataset('wall/2D/z', data=w['z'])
