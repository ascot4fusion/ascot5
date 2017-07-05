from pylab import *
import numpy as np

# The input file contains an integer giving the number of wall triangles
# followed by a list of coordinates for the corners of each triangle
# (x1,y1,z1,x2,y2,z2,x3,y3,z3).
def read_wall_3d(fn):
    
    f = open(fn,'r')

    str = {'n' : f.readline()}
    data = loadtxt(f)
    f.close()

    str['x1x2x3'] = [data[:,0], data[:,3], data[:,6]]
    str['y1y2y3'] = [data[:,1], data[:,4], data[:,7]]
    str['z1z2z3'] = [data[:,2], data[:,5], data[:,8]]

    return str


def read_wall_3d_ascot4(fn):
    
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
    f.readline()

    # Read data for each sector
    n_read = 0
    
    str['x1x2x3'] = np.array(np.zeros(num_lines, 3))
    str['y1y2y3'] = np.array(np.zeros(num_lines, 3))
    str['z1z2z3'] = np.array(np.zeros(num_lines, 3))
    for sector in range(str['n_sectors']):
        n_elements = int(f.readline().split()[0])
        f.readline()
        for i in range(n_elements):
            f.readline() # Skip empty line
            f.readline() # Skip element info
            for j in range(3):
                line = f.readline().split()
                str['x1x2x3'][n_read,j] = line[0]
                str['y1y2y3'][n_read,j] = line[1]
                str['z1z2z3'][n_read,j] = line[2]
                str['id'][n_read,j] = sector
            n_read = n_read + 1
            
    # fprintf('Printing the types: ')
    # for t=1:w.ntype
    #     fprintf('-');
    #     fprintf(fid,'\n%d elements belonging to sector type \n',w.types{t}.nelement);
    #     fprintf(fid,'%d  - number of vertices, divertor flag, ID, and vertices (x,y,z) \n',t);
    #     for e=1:w.types{t}.nelement
    #         fprintf(fid,'\n%d %d %s\n',...
    #             w.types{t}.vertexnum(e),...
    #             w.types{t}.divflag(e),...
    #             w.types{t}.id{e}          );
    #         % Here, the coordinates have to be written with sufficient floating
    #         % point accuracy. Otherwise it may turn out that two vertex points
    #         % appear to have the same coordinates which causes problems.
    #             fprintf(fid,'%5.7f %5.7f %5.7f\n', w.types{t}.coords{e}');
    #     end
    #     end
    f.close()
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
