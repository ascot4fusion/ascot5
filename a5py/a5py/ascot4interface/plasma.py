import numpy as np

def read_plasma(fn):
    with open(fn,'r') as f:
        c1 = f.readline()
        c2 = f.readline()
        c3 = f.readline()
        l = f.readline().split()
        if(fn[-2:] == '1d'):
            f.seek(0)
            return read_1d(f)
        elif(fn[-2:] == '2d'):
            f.seek(0)
            return read_2d(f)
        else:
            print('unrecognized first line after comments')
            print(l)

def read_1d(fh):
    pls = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nrho,nion = [int(number) for number in fh.readline().split()[:2]]
    pls['znum'] = np.array([int(znum) for znum in fh.readline().split()[:nion]])
    pls['anum'] = np.array([int(anum) for anum in fh.readline().split()[:nion]])
    pls['coll'] = np.array([int(coll) for coll in fh.readline().split()[:nion+1]])
    pls['nrho'] = nrho
    pls['nion'] = nion
    fh.readline() # ignore headers
    h5data = np.loadtxt(fh)
    pls['rho'] = np.array(h5data[:,0])
    pls['te'] = np.array(h5data[:,1])
    pls['ne'] = np.array(h5data[:,2])
    pls['vtor'] = np.array(h5data[:,3])
    pls['ti1'] = np.array(h5data[:,4])
    for i in range(0,nion):
        pls['ni'+str(i+1)] = np.array(h5data[:,5+i])
    return pls

def read_2d(fh):
    data = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nr,nz = [int(number) for number in fh.readline().split()[:2]]
    rmin,rmax,zmin,zmax = [float(number) for number in fh.readline().split()[:4]]
    nion,rho2d = [float(number) for number in fh.readline().split()[:2]]
    nion = int(nion)
    data['znum'] = np.array([float(znum) for znum in fh.readline().split()[:nion]])
    data['anum'] = np.array([float(anum) for anum in fh.readline().split()[:nion]])
    data['coll'] = np.array([float(coll) for coll in fh.readline().split()[:nion+1]])
    fieldnames = fh.readline().split()[0:-1:2]
    h5data = loadtxt(fh)
    data['r'] = linspace(rmin,rmax,nr)
    data['z'] = linspace(zmin,zmax,nz)
    for i,name in enumerate(fieldnames):
        data[name.lower()] = h5data[:,i].reshape(nr,nz).T

    if(i+1 < shape(h5data)[1]):
        print('Not all data fields assigned to struct')
        print('Something is probably wrong with the data file')

    return data

def write_plasma_1d(fn, p):
    dens_i = np.array([p['ni'+str(i)] for i in range(1,p['nion']+1)])
    p1D.write_hdf5(fn, p['nrho'], p['nion'], p['znum'], p['znum'], p['rho'],
                   np.zeros(p['rho'].shape), np.zeros(p['rho'].shape),
                   p['ne'], p['te'], dens_i, p['ti1'])
