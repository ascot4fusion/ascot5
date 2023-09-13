import h5py
import numpy as np

# The input file contains an integer giving the number of wall triangles
# followed by a list of coordinates for the corners of each triangle
# (x1,y1,z1,x2,y2,z2,x3,y3,z3).

def read_wall_3d(fn):

    # Crudely check number of lines to get maximum size for data
    num_lines = sum(1 for line in open(fn))
    copysector = False

    with open(fn,'r') as f:

        # fprintf(fid,'%d wall sectors; number, type, min. & max. phi(deg): \n\n',w.nsector);
        data = {}
        data['n_sectors'] = int(f.readline().split()[0])

        # Read sector ids
        data['ids'] = {}
        f.readline()
        for sector in range(data['n_sectors']):
            line = f.readline().split()
            data['ids'][int(float(line[0]))] = int(float(line[1]))

        # Read data for each sector
        n_read = 0

        data['x1x2x3'] = np.array(np.zeros((num_lines, 3)))
        data['y1y2y3'] = np.array(np.zeros((num_lines, 3)))
        data['z1z2z3'] = np.array(np.zeros((num_lines, 3)))
        data['id'] = np.array(np.zeros((num_lines, 1)))
        data['flag'] = np.array(np.zeros((num_lines, 1))) # This is not filled (for lack of test-files)
        for sector in range(1,data['n_sectors']+1):
            f.readline() # Skip empty line
            try:
                n_elements = int(float(f.readline().split()[0]))
            except Exception:
                print("Warning: There are multiple sectors but\n"
                      + "data only for one sector. The sector is\n"
                      + "copied and rotated." )
                copysector = True
                break
            print(f.readline())
            for i in range(n_elements):
                f.readline() # Skip empty line
                f.readline() # Skip element info
                for j in range(3):
                    line = f.readline().split()
                    data['x1x2x3'][n_read,j] = line[0]
                    data['y1y2y3'][n_read,j] = line[1]
                    data['z1z2z3'][n_read,j] = line[2]
                data['id'][n_read] = data['ids'][sector]
                n_read = n_read + 1

                # Some data may contain rectangles instead of triangles.
                # This just means we read the "extra" vertice and make another
                # triangle from it.
                pos = f.tell()
                line = f.readline()
                if len(line) == 0 or len(line) == 1:
                    # This was an empty line (no vertice data). Return to
                    # previous line.
                    f.seek(pos)
                else:
                    line = line.split()
                    data['x1x2x3'][n_read,0] = line[0]
                    data['y1y2y3'][n_read,0] = line[1]
                    data['z1z2z3'][n_read,0] = line[2]
                    data['x1x2x3'][n_read,1] = data['x1x2x3'][n_read-1,0]
                    data['y1y2y3'][n_read,1] = data['y1y2y3'][n_read-1,0]
                    data['z1z2z3'][n_read,1] = data['z1z2z3'][n_read-1,0]
                    data['x1x2x3'][n_read,2] = data['x1x2x3'][n_read-1,2]
                    data['y1y2y3'][n_read,2] = data['y1y2y3'][n_read-1,2]
                    data['z1z2z3'][n_read,2] = data['z1z2z3'][n_read-1,2]
                    data['id'][n_read] = data['ids'][sector]
                    n_read = n_read + 1

    # Remove extra zeros
    data['x1x2x3'] = data['x1x2x3'][0:n_read,:]
    data['y1y2y3'] = data['y1y2y3'][0:n_read,:]
    data['z1z2z3'] = data['z1z2z3'][0:n_read,:]
    data['id'] = data['id'][0:n_read,:]

    if copysector:
        r1r2r3 = np.sqrt(data['x1x2x3']*data['x1x2x3'] + data['y1y2y3']*data['y1y2y3'])
        t1t2t3 = np.arctan2(data['y1y2y3'], data['x1x2x3'])
        z1z2z3 = data['z1z2z3']
        data['x1x2x3'] = np.zeros( (n_read*data['n_sectors'],3) )
        data['y1y2y3'] = np.zeros( (n_read*data['n_sectors'],3) )
        data['z1z2z3'] = np.zeros( (n_read*data['n_sectors'],3) )

        for i in range(data['n_sectors']):
            x1x2x3 = r1r2r3 * np.cos(t1t2t3 + i*2*np.pi/data['n_sectors'])
            y1y2y3 = r1r2r3 * np.sin(t1t2t3 + i*2*np.pi/data['n_sectors'])
            data['x1x2x3'][n_read*i:n_read*(i+1),:] = x1x2x3
            data['y1y2y3'][n_read*i:n_read*(i+1),:] = y1y2y3
            data['z1z2z3'][n_read*i:n_read*(i+1),:] = z1z2z3

        data['id'] = np.linspace(1, n_read*data['n_sectors'],n_read*data['n_sectors'] )

    return data

def read_wall_3d_hdf5(fname):

    data = dict()
    
    with h5py.File(fname, 'r') as f: # Open for reading
        data['x1x2x3'] = f['/wall/3d/triangles_x1x2x3'][:]
        data['y1y2y3'] = f['/wall/3d/triangles_y1y2y3'][:]
        data['z1z2z3'] = f['/wall/3d/triangles_z1z2z3'][:]
        data[ "flag" ] = f['/wall/3d/triangles_flag'  ][:]
        
    data['n']=[[ data['x1x2x3'].shape[0] ]]

    return data
