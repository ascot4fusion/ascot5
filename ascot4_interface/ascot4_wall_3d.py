from pylab import *
import numpy as np

# The input file contains an integer giving the number of wall triangles
# followed by a list of coordinates for the corners of each triangle
# (x1,y1,z1,x2,y2,z2,x3,y3,z3).
def read_wall_3d(fn):
    
    f = open(fn,'r')

    data = loadtxt(f,skiprows=1)
    f.close()

    str = {}
    str['x1x2x3'] = [data[:,0], data[:,3], data[:,6]]
    str['y1y2y3'] = [data[:,1], data[:,4], data[:,7]]
    str['z1z2z3'] = [data[:,2], data[:,5], data[:,8]]

    return str

def write_wall_3d(f, w):
    f.create_dataset('wall/3D/x1x2x3', data=w['x1x2x3'])
    f.create_dataset('wall/3D/y1y2y3', data=w['y1y2y3'])
    f.create_dataset('wall/3D/z1z2z3', data=w['z1z2z3'])
    f['wall/3D'].attrs['min_x'] = np.amin(w['x1x2x3'])
    f['wall/3D'].attrs['max_x'] = np.amax(w['x1x2x3'])
    f['wall/3D'].attrs['min_y'] = np.amin(w['y1y2y3'])
    f['wall/3D'].attrs['max_y'] = np.amax(w['y1y2y3'])
    f['wall/3D'].attrs['min_z'] = np.amin(w['z1z2z3'])
    f['wall/3D'].attrs['max_z'] = np.amax(w['z1z2z3'])
