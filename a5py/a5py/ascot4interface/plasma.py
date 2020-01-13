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
    data = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nrho,nion = [int(float(number)) for number in fh.readline().split()[:2]]
    data['znum'] = np.array([int(float(znum)) for znum in fh.readline().split()[:nion]])
    data['anum'] = np.array([int(float(anum)) for anum in fh.readline().split()[:nion]])
    data['coll'] = np.array([int(float(coll)) for coll in fh.readline().split()[:nion+1]])
    data['nrho'] = nrho
    data['nion'] = nion
    fh.readline() # ignore headers
    h5data = np.loadtxt(fh)
    data['rho'] = np.array(h5data[:,0])
    data['te'] = np.array(h5data[:,1])
    data['ne'] = np.array(h5data[:,2])
    data['vtor'] = np.array(h5data[:,3])
    data['ti'] = np.array(h5data[:,4])
    for i in range(0,nion):
        data['ni'+str(i+1)] = np.array(h5data[:,5+i])
    # Make sure the input is linearly spaced. If not, interpolate
    tol = 1.0001
    diff = np.diff(data['rho'])
    if ( max(diff)/min(diff) > tol):
        print("Warning! Interpolating plasma data to uniform grid")
        new_rho = np.linspace(np.amin(data['rho']),
                              np.amax(data['rho']),
                              data['nrho'])
        data['ne'] = np.interp(new_rho, data['rho'], data['ne'])
        data['te'] = np.interp(new_rho, data['rho'], data['te'])
        for i in range(1, data['nion']+1):
            data['ni'+str(i)] = np.interp(new_rho, data['rho'],
                                          data['ni'+str(i)])
        data['ti'] = np.interp(new_rho, data['rho'], data['ti'])
        data['rho'] = new_rho
    data['ni'] = np.array([data['ni'+str(i)] for i in range(1,data['nion']+1)])
    data['ni'] = np.transpose(data['ni'])
    # Add extra data point outside rho=1 to avoid out of data range errors
    if ( np.amax(data['rho']) <= 1.0 ):
        print("Warning! Adding small datapoint outside rho=1.0")
        data['nrho'] = data['nrho'] + 1
        data['rho'] = np.append(data['rho'], 2*data['rho'][-1]-data['rho'][-2])
        data['ne'] = np.append(data['ne'], data['ne'][-1]*1e-10)
        data['te'] = np.append(data['te'], data['te'][-1])
        data['ni'] = np.append(data['ni'],
                               np.expand_dims(data['ni'][-1,:]*1e-10, 1).T,
                               axis=0)
        data['ti'] = np.append(data['ti'], data['ti'][-1])
    return data

def read_2d(fh):
    data = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nr,nz = [int(float(number)) for number in fh.readline().split()[:2]]
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
                   p['ne'], p['te'], dens_i, p['ti'])
