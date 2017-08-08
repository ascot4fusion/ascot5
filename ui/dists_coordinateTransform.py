import numpy as np
#import scipy.interpolate as sp
from scipy.interpolate import griddata

def vpavpe2Epitch(a5dist, E_edges, xi_edges, mass):
    
    # Check that vpa and vpe are present and note what 
    # indices they correspond to
    Ivpa = a5dist['ordinate_names'].index("vpa")
    Ivpe = a5dist['ordinate_names'].index("vpe")

    # Transform given E-xi grid to vpa vpe coordinates
    E  = np.linspace(0.5*(E_edges[0]+E_edges[1]), 0.5*(E_edges[-2]+E_edges[-1]), num=E_edges.size-1)
    xi = np.linspace(0.5*(xi_edges[0]+xi_edges[1]), 0.5*(xi_edges[-2]+xi_edges[-1]), num=xi_edges.size-1)

    v = np.sqrt( 2 * E * 1.6022e-19 / (mass * 1.660539040e-27) )
    
    vg, xig = np.meshgrid(v,xi)

    vpag = xig * vg
    vpeg = np.sqrt( 1 - xig * xig ) * vg

    points = np.transpose(np.array([vpag.flatten(), vpeg.flatten()]))

    vpagrid, vpegrid = np.meshgrid(a5dist['vpa'], a5dist['vpe'])
    grid = np.transpose(np.array([vpagrid.flatten(), vpegrid.flatten()]))

    vg = vg.flatten()
    xig = xig.flatten()

    # Interpolate distribution at desired coordinates
    f = np.zeros((a5dist['nR'], a5dist['nz'], E.size, xi.size, a5dist['ntime']))
    for iR in range(0, a5dist['nR']):
        for iz in range(0, a5dist['nz']):
            for itime in range(0, a5dist['ntime']):
                a = griddata(grid, np.squeeze(a5dist['ordinate'][iR,iz,:,:,itime]).flatten(), points, method='linear', fill_value=0)
                a = a * vg / np.sqrt(1 - np.power(xig,2))
                a = (np.reshape(a,(xi.size,E.size),'C'))
                f[iR,iz,:,:,itime] = np.transpose( a ) #* ((a5dist['vpa'][1]-a5dist['vpa'][0])*(a5dist['vpe'][1]-a5dist['vpe'][0])) / ((v[1]-v[0])*(xi[1]-xi[0]))


    # Create E-xi distribution
    Exidist = {}
    Exidist['ordinate'] = f
    Exidist['E'] = E
    Exidist['xi'] = xi

    return Exidist
