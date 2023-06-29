import numpy as np

def read_neutral(fn):
    with open(fn,'r') as f:
        if(fn[-2:] == '1d'):
            f.seek(0)
            return read_1d(f)
        elif(fn[-2:] == '3d'):
            print("3D neutral data has not yet been implemented for ASCOT4, "
                  "and neither has a method for reading it here.")
        else:
            print('Unrecognized suffix in filename')

def read_1d(fh):
    data = {'comm1' : fh.readline(),
            'comm2' : fh.readline(),
            'comm3' : fh.readline()}
    nrho = int(float(fh.readline()))
    nspecies = 1 # ASCOT4 does not support several neutral species
    data['nrho'] = nrho
    data['nspecies'] = nspecies
    fh.readline() # ignore headers
    h5data = np.loadtxt(fh)
    data['rho'] = np.array(h5data[:,0])
    for i in range(0, data['nspecies']):
        data['dens'+str(i+1)] = np.array(h5data[:,1+i])
    data['temp'] = np.array(h5data[:,2])
    # Make sure the input is linearly spaced. If not, interpolate
    tol = 1.0001
    diff = np.diff(data['rho'])
    if ( max(diff)/min(diff) > tol):
        print("Warning! Interpolating neutral data to uniform grid")
        new_rho = np.linspace(np.amin(data['rho']),
                              np.amax(data['rho']),
                              data['nrho'])
        for i in range(0, data['nspecies']):
            data['dens'+str(i+1)] = np.interp(new_rho, data['rho'],
                                              data['dens'+str(i+1)])
        data['temp'] = np.interp(new_rho, data['rho'], data['temp'])
        data['rho'] = new_rho
    data['dens'] = np.array([data['dens'+str(i+1)] for i in range(0, data['nspecies'])])
    data['dens'] = np.transpose(data['dens'])
    # Add extra data point outside rho=1 to avoid out of data range errors
    if ( np.amax(data['rho']) <= 1.0 ):
        print("Warning! Adding small datapoint outside rho=1.0")
        data['nrho'] = data['nrho'] + 1
        data['rho'] = np.append(data['rho'], 2*data['rho'][-1]-data['rho'][-2])
        data['dens'] = np.append(data['dens'],
                                 np.expand_dims(data['dens'][-1,:]*1e-10, 1).T,
                                 axis=0)
        data['temp'] = np.append(data['temp'], data['temp'][-1])
    return data
