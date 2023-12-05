import numpy as np
import unyt
from scipy.interpolate import griddata
from scipy.io.netcdf import netcdf_file as Dataset
import h5py
import matplotlib.pyplot as plt

from a5py.ascot5io.dist import DistData

def transp_exidist(transpCDF, ascotH5, nr=30, nz=42, plot=True):
    """Convert TRANSP distribution to ASCOT5 format.
    """
    # DISTRIBUTION
    distName = 'rzPitchEdist'
    tra = transp2dict(transpCDF)

    nE = len(tra['energy'])
    nP = len(tra['pitch'])
    nT = 1
    nSpecies = 1
    nOrdDim = 1

    def gridminmax(x):
        xmin = np.amin(x)
        xmax = np.amax(x)
        return xmin - 0.1 * ( xmax - xmin ), xmax + 0.1 * ( xmax - xmin )

    rmin, rmax = gridminmax(tra["grid"]["r"])
    zmin, zmax = gridminmax(tra["grid"]["z"])

    dist = DistData(np.zeros((nr, nz, nP, nE)),
                    r=np.linspace(rmin, rmax, nr + 1) * unyt.m,
                    z=np.linspace(zmin, zmax, nz + 1) * unyt.m,
                    pitch=tra["pitch_boundary"]  * unyt.dimensionless,
                    ekin=tra["energy_boundary"] * unyt.eV)

    R, Z = np.meshgrid(dist.abscissa("r"), dist.abscissa("z"))
    for ip in range(nP):
        for ie in range(nE):
            dist._distribution[:,:, ip, ie] = griddata(
                (tra["grid"]["r"], tra["grid"]["z"]),
                tra["dist"][:,ip,ie], (R, Z), method='linear', fill_value=0.0).T
            
    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(2,2,1)
        ax2 = fig.add_subplot(2,2,2)
        dist.integrate(copy=True, pitch=np.s_[:], ekin=np.s_[:]).plot(axes=ax1)
        dist.integrate(copy=True, r=np.s_[:], z=np.s_[:]).plot(axes=ax2)

        ax3 = fig.add_subplot(2,2,3)
        ax4 = fig.add_subplot(2,2,4)
        h = ax3.scatter(tra['grid']['r'], tra['grid']['z'], s=50, c=tra['n_fi'],
                    edgecolors='k', linewidths=0.5)
        plt.colorbar(h,ax=ax3)
        h = ax4.pcolormesh(tra['pitch_boundary'], tra['energy_boundary'],
                           np.sum(tra['dist'], axis=0).T)
        plt.colorbar(h,ax=ax4)

        plt.show()

def transp2dict(filename, plot=False):
    tra = {}
    
    with Dataset(filename, "r", mmap=False) as f:
        ncfile = f.variables.copy()
        tra['filename'] = filename
        tra['grid'] = {
            'z': ncfile['Z2D'][:].copy() / 100,
            'r': ncfile['R2D'][:].copy() / 100,
            'theta': ncfile['TH2D'][:].copy(),
            'X': ncfile['X2D'][:].copy(),
            'volumes': ncfile['BMVOL'][:].copy()
        }

        tra['pitch']           = ncfile['A_D_NBI'][:].copy()
        tra['pitch_boundary']  = ncfile['AB_D_NBI'][:].copy()
        tra['energy']          = ncfile['E_D_NBI'][:].copy()
        tra['energy_boundary'] = ncfile['EB_D_NBI'][:].copy()
        tra['energy_units']    = 'eV'
        tra['dist']            = ncfile['F_D_NBI'][:].copy() * 0.5
        tra['dist_units']      = ncfile['F_D_NBI'][:].copy()
        tra['time']            = ncfile['TIME'][()].copy()
        tra['runid']           = ncfile['TRANSP_RUNID'][:].copy()

    tra['n_fi'] = np.squeeze(np.sum(np.sum(tra['dist'], axis=2), axis=1))
    tra['n_fi'] *= 100**3  # 1/cm^3 => 1/m^3
    tra['n_fi'] *= (tra['pitch_boundary'][1] - tra['pitch_boundary'][0])
    tra['n_fi'] *= (tra['energy_boundary'][1] - tra['energy_boundary'][0])
    tra['volumeAverage'] = np.sum( tra['n_fi'] * tra['grid']['volumes'] ) \
        / np.sum(tra['grid']['volumes'])
    vol = np.transpose(np.tile(tra['grid']['volumes'], (75, 50, 1)))
    tra['vdist'] = np.sum(tra['dist'] * vol, axis=2)

    if plot:
        r = np.linspace(np.min(tra['grid']['r']), np.max(tra['grid']['r']), 30)
        z = np.linspace(np.min(tra['grid']['z']), np.max(tra['grid']['z']), 50)

        d = griddata((tra['grid']['r'], tra['grid']['z']), tra['n_fi'],
                 (r[None, :], z[:, None]), method='linear')
    
        plt.pcolormesh(r, z, d)
        plt.scatter(tra['grid']['r'], tra['grid']['z'], s=50, c=tra['n_fi'],
                    edgecolors='k', linewidths=0.5)
        plt.axis('image')
        plt.xlabel('R (m)')
        plt.ylabel('z (m)')
        cbh = plt.colorbar()
        cbh.set_label('1/(m^3)')
        plt.title('Fast ion density')
        plt.show()

    return tra

def convert_ascot4_4dDist_PE_2_VpaVpe(a):
    # Replace with actual implementation
    pass

def write_ascot4_distribution(h5file, data, distName):
    # Replace with actual implementation
    pass

def write_ascot4_species(h5file, species):
    # Replace with actual implementation
    pass