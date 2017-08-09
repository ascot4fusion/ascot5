import h5py
from pylab import *
import numpy as np

# # The input file contains an integer giving the number of wall triangles
# # followed by a list of coordinates for the corners of each triangle
# # (x1,y1,z1,x2,y2,z2,x3,y3,z3).
# def read_wall_3d(fn):
    
#     f = open(fn,'r')

#     str = {'n' : f.readline()}
#     data = loadtxt(f)
#     f.close()

#     str['x1x2x3'] = [data[:,0], data[:,3], data[:,6]]
#     str['y1y2y3'] = [data[:,1], data[:,4], data[:,7]]
#     str['z1z2z3'] = [data[:,2], data[:,5], data[:,8]]

#     return str

def read_wall_3d(fn):
    
    # Crudely check number of lines to get maximum size for data
    num_lines = sum(1 for line in open(fn))

    f = open(fn,'r')

    # fprintf(fid,'%d wall sectors; number, type, min. & max. phi(deg): \n\n',w.nsector);
    str = {}
    str['n_sectors'] = int(f.readline().split()[0])

    # Read sector ids
    str['ids'] = {}
    f.readline()
    for sector in range(str['n_sectors']):
        line = f.readline().split()
        str['ids'][int(line[0])] = int(line[1])

    # Read data for each sector
    n_read = 0
    
    str['x1x2x3'] = np.array(np.zeros((num_lines, 3)))
    str['y1y2y3'] = np.array(np.zeros((num_lines, 3)))
    str['z1z2z3'] = np.array(np.zeros((num_lines, 3)))
    str['id'] = np.array(np.zeros((num_lines, 1)))
    for sector in range(1,str['n_sectors']+1):
        f.readline() # Skip empty line
        n_elements = int(f.readline().split()[0])
        print f.readline()
        for i in range(n_elements):
            f.readline() # Skip empty line
            f.readline() # Skip element info
            for j in range(3):
                line = f.readline().split()
                str['x1x2x3'][n_read,j] = line[0]
                str['y1y2y3'][n_read,j] = line[1]
                str['z1z2z3'][n_read,j] = line[2]
            str['id'][n_read] = str['ids'][sector]
            n_read = n_read + 1
            
    f.close()

    # Remove extra zeros
    str['x1x2x3'] = str['x1x2x3'][0:n_read,:]
    str['y1y2y3'] = str['y1y2y3'][0:n_read,:]
    str['z1z2z3'] = str['z1z2z3'][0:n_read,:]
    str['id'] = str['id'][0:n_read,:]    
    
    return str

def read_wall_3d_hdf5(f):
    
    # fprintf(fid,'%d wall sectors; number, type, min. & max. phi(deg): \n\n',w.nsector);
    str = {}
    str['n_sectors'] = int(f.readline().split()[0])

    # Read sector ids
    str['ids'] = {}
    f.readline()
    for sector in range(str['n_sectors']):
        line = f.readline().split()
        str['ids'][int(line[0])] = int(line[1])

    # Read data for each sector
    n_read = 0
    
    str['x1x2x3'] = np.array(np.zeros((num_lines, 3)))
    str['y1y2y3'] = np.array(np.zeros((num_lines, 3)))
    str['z1z2z3'] = np.array(np.zeros((num_lines, 3)))
    str['id'] = np.array(np.zeros((num_lines, 1)))
    for sector in range(1,str['n_sectors']+1):
        f.readline() # Skip empty line
        n_elements = int(f.readline().split()[0])
        print f.readline()
        for i in range(n_elements):
            f.readline() # Skip empty line
            f.readline() # Skip element info
            for j in range(3):
                line = f.readline().split()
                str['x1x2x3'][n_read,j] = line[0]
                str['y1y2y3'][n_read,j] = line[1]
                str['z1z2z3'][n_read,j] = line[2]
            str['id'][n_read] = str['ids'][sector]
            n_read = n_read + 1
            
    f.close()

    # Remove extra zeros
    str['x1x2x3'] = str['x1x2x3'][0:n_read,:]
    str['y1y2y3'] = str['y1y2y3'][0:n_read,:]
    str['z1z2z3'] = str['z1z2z3'][0:n_read,:]
    str['id'] = str['id'][0:n_read,:]    
    
    f.close()

    return str

def write_wall_3d(f, w):
    f.create_dataset('wall/3D/x1x2x3', data=w['x1x2x3'])
    f.create_dataset('wall/3D/y1y2y3', data=w['y1y2y3'])
    f.create_dataset('wall/3D/z1z2z3', data=w['z1z2z3'])
    f.create_dataset('wall/3D/flag', data=w['id'])
    f['wall/3D'].attrs['n_elements'] = size(w['id'])
    f['wall/3D'].attrs['min_x'] = np.amin(w['x1x2x3'])
    f['wall/3D'].attrs['max_x'] = np.amax(w['x1x2x3'])
    f['wall/3D'].attrs['min_y'] = np.amin(w['y1y2y3'])
    f['wall/3D'].attrs['max_y'] = np.amax(w['y1y2y3'])
    f['wall/3D'].attrs['min_z'] = np.amin(w['z1z2z3'])
    f['wall/3D'].attrs['max_z'] = np.amax(w['z1z2z3'])

def write_wall_3d_hdf5(f, inputf):
    f.create_group('wall/3D')
    f.create_dataset('wall/3D/x1x2x3', data=inputf["/wall/3d/triangles_x1x2x3"])
    f.create_dataset('wall/3D/y1y2y3', data=inputf["/wall/3d/triangles_y1y2y3"])
    f.create_dataset('wall/3D/z1z2z3', data=inputf["/wall/3d/triangles_z1z2z3"])
    f.create_dataset('wall/3D/flag', data=inputf["/wall/3d/triangles_flag"])
    f['wall/3D'].attrs['n_elements'] = size(inputf["/wall/3d/triangles_flag"])
    f['wall/3D'].attrs['min_x'] = np.amin(inputf["/wall/3d/triangles_x1x2x3"])
    f['wall/3D'].attrs['max_x'] = np.amax(inputf["/wall/3d/triangles_x1x2x3"])
    f['wall/3D'].attrs['min_y'] = np.amin(inputf["/wall/3d/triangles_y1y2y3"])
    f['wall/3D'].attrs['max_y'] = np.amax(inputf["/wall/3d/triangles_y1y2y3"])
    f['wall/3D'].attrs['min_z'] = np.amin(inputf["/wall/3d/triangles_z1z2z3"])
    f['wall/3D'].attrs['max_z'] = np.amax(inputf["/wall/3d/triangles_z1z2z3"])
