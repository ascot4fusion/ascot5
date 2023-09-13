import numpy as np

def read_erad(fn):

    with open(fn,'r') as f:
        data = {'comm1' : f.readline()} #Skip comment line
        data['n_rho'] = int(float(f.readline().split()[0]))

        h5data = np.loadtxt(f)

        data['rho'] = h5data[:,0]
        data['dV_drho'] = h5data[:,1]
        # For data in format dV/rho, we can ignore effective minor radius

    # Make sure the input is linearly spaced. If not, interpolate
    tol = 1.0001
    diff = np.diff(data['rho'])
    if ( max(diff)/min(diff) > tol):
        print("Warning! Interpolating dV_drho to uniform grid")
        new_rho = np.linspace(np.amin(data['rho']),
                              np.amax(data['rho']),
                              data['n_rho'])
        data['dV_drho'] = np.interp(new_rho, data['rho'],
                                    data['dV_drho'])
        data['rho'] = new_rho

    return data

def write_ascot4_input_erad(filename,plasmaData,comment):
    '''
    Write the radial electric field to input.erad
    @param plasmaData: dictionary with keys:
        'rho','Er'
    '''

    header=comment+ '\n'

    Nrad = len(plasmaData['rho'])
    header = header + ' {} # Nradial points\n'.format(Nrad)


    D=np.vstack( (plasmaData['rho'],plasmaData['Er']) )
    np.savetxt(filename,np.transpose(D),delimiter=' ',comments='',header=header)
