import numpy as np

def read(dists, mode):
    out = {}

    if "rzVDist" in dists:
        out['R_edges']    = dists['rzVDist/abscissae/dim1'][:]
        out['R']          = np.linspace(0.5*(out['R_edges'][0]+out['R_edges'][1]), 0.5*(out['R_edges'][-2]+out['R_edges'][-1]), num=out['R_edges'].size-1)
        out['R_unit']     = 'm'
        out['nR']         = out['R'].size

        out['z_edges']    = dists['rzVDist/abscissae/dim2'][:]
        out['z']          = np.linspace(0.5*(out['z_edges'][0]+out['z_edges'][1]), 0.5*(out['z_edges'][-2]+out['z_edges'][-1]), num=out['z_edges'].size-1)
        out['z_unit']     = 'm'
        out['nz']         = out['z'].size

        out['vpa_edges']  = dists['rzVDist/abscissae/dim3'][:]
        out['vpa']        = np.linspace(0.5*(out['vpa_edges'][0]+out['vpa_edges'][1]), 0.5*(out['vpa_edges'][-2]+out['vpa_edges'][-1]), num=out['vpa_edges'].size-1)
        out['vpa_unit']   = 'm/s'
        out['nvpa']       = out['vpa'].size

        out['vpe_edges']  = dists['rzVDist/abscissae/dim4'][:]
        out['vpe']        = np.linspace(0.5*(out['vpe_edges'][0]+out['vpe_edges'][1]), 0.5*(out['vpe_edges'][-2]+out['vpe_edges'][-1]), num=out['vpe_edges'].size-1)
        out['vpe_unit']   = 'm/s'
        out['nvpe']       = out['vpe'].size

        out['time_edges'] = dists['rzVDist/abscissae/dim5'][:]
        out['time']       = np.linspace(0.5*(out['time_edges'][0]+out['time_edges'][1]), 0.5*(out['time_edges'][-2]+out['time_edges'][-1]), num=out['time_edges'].size-1)
        out['time_unit']  = 's'
        out['ntime']      = out['time'].size

        out['ordinate']       = np.reshape(dists['rzVDist/ordinate'][:].flatten(), [out['nR'], out['nz'], out['nvpa'], out['nvpe'], out['ntime']] )
        out['ordinate_names'] = ['R','z','vpa','vpe','time']

    return out
