"""
Time-dependent Non-axisymmetric tokamak magnetic field HDF5 IO

File: B_3DST.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
               b_phimin, b_phimax, b_nphi, b_tmin, b_tmax, b_nt,
               axisr, axisz, psi, psi0, psi1, br, bphi, bz,
               psi_rmin=None, psi_rmax=None, psi_nr=None,
               psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
    """
    Write 3DST magnetic field input in HDF5 file.

    Note that br and bz should not include the equilibrium component of the
    magnetic field as that is calculated from psi by ASCOT5 during the
    simulation.

    It is possible to use different Rz grids for psi and magnetic field
    components by giving Rz grid for psi separately.

    The toroidal angle phi is treated as a periodic coordinate meaning that
    B(phi) = B(phi + n*(b_phimax - b_phimin)). Do note that to avoid duplicate
    data, the last points in phi axis in B data are not at b_phimax, i.e.
    br[:,:,-1] != BR(phi=b_phimax).

    Args:
        fn : str <br>                                                               
            Full path to the HDF5 file.                                             
        b_rmin : float <br>                                                         
            Magnetic field data R grid min edge [m].                                
        b_rmax : float <br>                                                         
            Magnetic field data R grid max edge [m].                                
        b_nr : int <br>                                                             
            Number of R grid points in magnetic field data.                         
        b_zmin : float <br>                                                         
            Magnetic field data z grid min edge [m].                                
        b_zmax : float <br>                                                         
            Magnetic field data z grid max edge [m].                                
        b_nz : int <br>                                                             
            Number of z grid points in magnetic field data.                         
        b_phimin : float <br>                                                       
            Magnetic field data phi grid min edge [deg].                            
        b_phimax : float <br>                                                       
            Magnetic field data phi grid max edge [deg].                            
        b_nphi : int <br>                                                           
            Number of phi grid points in magnetic field data.                       
        b_tmin : float <br>                                                       
            Magnetic field data time grid min edge [s].                            
        b_tmax : float <br>                                                       
            Magnetic field data time grid max edge [s].                            
        b_nt : int <br>                                                           
            Number of t grid points in magnetic field data.                       
        axisr : float <br>                                                          
            Magnetic axis R coordinate [m].                                         
        axisz : float <br>                                                          
            Magnetic axis z coordinate [m].                                         
        psi0 : float <br>                                                           
            On-axis poloidal flux value [Vs/m].                                     
        psi1 : float <br>                                                           
            Separatrix poloidal flux value [Vs/m].                                  
        psi : array_like (nr, nz) <br>                                              
            Poloidal flux values on the Rz grid [Vs/m].                             
        br : array_like (nr,nphi,nz) <br>                                           
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
            Matrices shape must be (n_time,n_z,n_phi,n_R) CHECK THIS
        bphi : array_like (nr,nphi,nz) <br>                                         
            Magnetic field phi component on Rz grid [T].
            Matrices shape must be (n_time,n_z,n_phi,n_R) CHECK THIS
        bz : array_like (nr,nphi,nz) <br>                                           
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
            Matrices shape must be (n_time,n_z,n_phi,n_R) CHECK THIS
        psi_rmin : float, optional <br>                                             
            Psi data R grid min edge [m].                                           
        psi_rmax : float, optional <br>                                             
            Psi data R grid max edge [m].                                           
        psi_nr : int, optional <br>                                                 
            Number of R grid points in psi data.                                    
        psi_zmin : float, optional <br>                                             
            Psi data z grid min edge [m].                                           
        psi_zmax : float, optional <br>                                             
            Psi data z grid max edge [m].                                           
        psi_nz : int, optional <br>                                                 
            Number of z grid points in psi data.                                    
        desc : str, optional <br>                                                   
            Input description.                

    Returns:
        Name of the new input that was written.

    Notes
    -------

    Rphiz coordinates form a right-handed cylindrical coordinates.
    phi is in degrees and is considered as a periodic coordinate.
    b_phimin is where the period begins and phimax is the last data point,
    i.e. not the end of period. E.g if you have a n = 2 symmetric system
    with nphi = 180 deg and b_phimin = 0 deg, then phimax should be 179 deg.

    Within ASCOT5, the magnetic field is evaluated as:

    br = br' + dPsi/dz,
    bphi = bphi',
    bz = bz' + dPsi/dR,

    where ' notates input fields.
    """

    parent = "bfield"
    group  = "B_3DST"
    gname  = ""                                                                     

    # Define psigrid to be same as Bgrid if not stated otherwise.
    if(psi_rmin is None or psi_rmax is None or psi_nr is None or
       psi_zmin is None or psi_zmax is None or psi_nz is None):
        psi_rmin = b_rmin
        psi_rmax = b_rmax
        psi_nr   = b_nr
        psi_zmin = b_zmin
        psi_zmax = b_zmax
        psi_nz   = b_nz

    ## continue here (and delete this comment)
        
    br = np.transpose(br,(0,2,1,3))
    bphi = np.transpose(bphi,(0,2,1,3))
    bz = np.transpose(bz,(0,2,1,3))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("R_min",         (1,), data=b_rmin,    dtype="f8")
        g.create_dataset("R_max",         (1,), data=b_rmax,    dtype="f8")
        g.create_dataset("n_R",           (1,), data=b_nr,      dtype="i8")
        g.create_dataset("phi_min",       (1,), data=b_phimin,  dtype="f8")
        g.create_dataset("phi_max",       (1,), data=b_phimax,  dtype="f8")
        g.create_dataset("n_phi",         (1,), data=b_nphi,    dtype="i8")
        g.create_dataset("z_min",         (1,), data=b_zmin,    dtype="f8")
        g.create_dataset("z_max",         (1,), data=b_zmax,    dtype="f8")
        g.create_dataset("n_z",           (1,), data=b_nz,      dtype="i8")
        g.create_dataset("t_min",         (1,), data=b_tmin,    dtype="f8")
        g.create_dataset("t_max",         (1,), data=b_tmax,    dtype="f8")
        g.create_dataset("n_t",           (1,), data=b_nt,   dtype="i8")
        g.create_dataset("psigrid_R_min", (1,), data=psi_rmin,   dtype="f8")
        g.create_dataset("psigrid_R_max", (1,), data=psi_rmax,   dtype="f8")
        g.create_dataset("psigrid_n_R",   (1,), data=psi_nr,     dtype="i8")
        g.create_dataset("psigrid_z_min", (1,), data=psi_zmin,   dtype="f8")
        g.create_dataset("psigrid_z_max", (1,), data=psi_zmax,   dtype="f8")
        g.create_dataset("psigrid_n_z",   (1,), data=psi_nz,     dtype="i8")
        g.create_dataset("psi",                 data=psi,   dtype="f8")
        g.create_dataset("br",                 data=br,     dtype="f8")
        g.create_dataset("bphi",               data=bphi,   dtype="f8")
        g.create_dataset("bz",                 data=bz,     dtype="f8")
        g.create_dataset("axis_R",        (1,), data=axisr,   dtype="f8")
        g.create_dataset("axis_z",        (1,), data=axisz,   dtype="f8")
        g.create_dataset("psi0",          (1,), data=psi0, dtype="f8")
        g.create_dataset("psi1",          (1,), data=psi1, dtype="f8")


def read_hdf5(fn, qid):
    """
    Read 3D magnetic field input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the bfield to be read.

    Returns
    -------

    Dictionary containing magnetic field data.
    """

    path = "bfield" + "/B_3DST-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["R_min"] = f[path]["R_min"][:]
        out["R_max"] = f[path]["R_max"][:]
        out["n_R"]   = f[path]["n_R"][:]

        out["phi_min"] = f[path]["phi_min"][:]
        out["phi_max"] = f[path]["phi_max"][:]
        out["n_phi"]   = f[path]["n_phi"][:]

        out["z_min"] = f[path]["z_min"][:]
        out["z_max"] = f[path]["z_max"][:]
        out["n_z"]   = f[path]["n_z"][:]

        out["t_min"] = f[path]["t_min"][:]
        out["t_max"] = f[path]["t_max"][:]
        out["n_t"]   = f[path]["n_t"][:]

        out["psi"]   = f[path]["psi"][:]
        out["br"]   = f[path]["br"][:]
        out["bphi"] = f[path]["bphi"][:]
        out["bz"]   = f[path]["bz"][:]

        out["axis_R"] = f[path]["axis_R"][:]
        out["axis_z"] = f[path]["axis_z"][:]

        out["psiaxis"] = f[path]["psi0"][:]
        out["psi1"] = f[path]["psi1"][:]

    return out

class B_3DST(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
