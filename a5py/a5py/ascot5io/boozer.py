"""
Boozer coordinate input IO.

File: boozer.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, psimin, psimax, npsi, ntheta, rmin, rmax, nr, zmin, zmax, nz,
               r0, z0, psi0, psi1, psi_rz, theta_psithetageom, nu_psitheta,
               nrzs, rs, zs, desc=None):
    """
    Write boozer input to HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        psimin : float <br>
            Minimum psi grid value.
        psimax : float <br>
            Maximum psi grid value.
        npsi : int <br>
            Number of psi grid points.
        ntheta : int <br>
            Number of theta grid (both boozer and geometric) values.
        rmin : float <br>
            Minimum R grid value.
        rmax : float <br>
            Maximum R grid value.
        nr : int <br>
            Number of R grid points.
        zmin : float <br>
            Minimum z grid value.
        zmax : float <br>
            Maximum z grid value.
        nz : int <br>
            Number of z grid points.
        r0 : float <br>
            Magnetic axis R coordinate.
        z0 : float <br>
            Magnetic axis z coordinate.
        psi0 : float <br>
            Coordinate psi on axis.
        psi1 : float <br>
            Coordinate psi on separatrix.
        psi_rz : array_like (nr,nz) <br>
            Coordinate psi(R,z).
        theta_psithetageom : array_like (npsi,ntheta) <br>
            Coordinate theta(psi, thetageom).
        nu_psitheta : array_like (npsi,ntheta) <br>
            nu(psi, theta).
        nrsz : int <br>
            Number of separatrix Rz points.
        rs : array_like (nrsz,1) <br>
            Separatrix R coordinates.
        zs : array_like (nrsz,1) <br>
            Separatrix z coordinates.
        desc : str, optional <br>
            Input's description.

    Returns:
        Name of the new input that was written.
    """
    assert psi_rz.shape             == (nr,nz)
    assert theta_psithetageom.shape == (npsi,ntheta)
    assert nu_psitheta.shape        == (npsi,ntheta)
    assert rs.size == nrzs
    assert zs.size == nrzs

    psi_rz             = np.transpose(psi_rz)
    theta_psithetageom = np.transpose(theta_psithetageom)
    nu_psitheta        = np.transpose(nu_psitheta)

    parent = "boozer"
    group  = "Boozer"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        # grid specifications
        g.create_dataset("psimin", (1,), data=psimin, dtype="f8")
        g.create_dataset("psimax", (1,), data=psimax, dtype="f8")
        g.create_dataset("npsi",   (1,), data=npsi,   dtype="i4")
        g.create_dataset("ntheta", (1,), data=ntheta, dtype="i4")
        g.create_dataset("rmin",   (1,), data=rmin,   dtype="f8")
        g.create_dataset("rmax",   (1,), data=rmax,   dtype="f8")
        g.create_dataset("nr",     (1,), data=nr,     dtype="i4")
        g.create_dataset("zmin",   (1,), data=zmin,   dtype="f8")
        g.create_dataset("zmax",   (1,), data=zmax,   dtype="f8")
        g.create_dataset("nz",     (1,), data=nz,     dtype="i4")
        g.create_dataset("r0",     (1,), data=r0,     dtype="f8")
        g.create_dataset("z0",     (1,), data=z0,     dtype="f8")
        g.create_dataset("nrzs",   (1,), data=nrzs,   dtype="i4")

        # the outermost poloidal psi-surface contour
        g.create_dataset("rs", (nrzs,), data=rs, dtype="f8")
        g.create_dataset("zs", (nrzs,), data=zs, dtype="f8")

        # psi data min and max values
        g.create_dataset("psi0", (1,), data=psi0, dtype="f8")
        g.create_dataset("psi1", (1,), data=psi1, dtype="f8")

        # tabulated coordinates maps
        g.create_dataset("psi_rz", (nz,nr), data=psi_rz, dtype="f8")
        g.create_dataset("theta_psithetageom", (ntheta,npsi),
                         data=theta_psithetageom, dtype="f8")
        g.create_dataset("nu_psitheta", (ntheta,npsi),
                         data=nu_psitheta, dtype="f8")

    return gname


def read_hdf5(fn, qid):
    """
    Read Boozer input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "boozer/Boozer_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    out["psi_rz"]             = np.transpose(out["psi_rz"])
    out["theta_psithetageom"] = np.transpose(out["theta_psithetageom"])
    out["nu_psitheta"]        = np.transpose(out["nu_psitheta"])
    return out


def write_hdf5_dummy(fn, desc="Dummy"):
    """
    Write dummy boozer input.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """

    psimin     = 0
    psimax     = 1
    npsi       = 6
    ntheta     = 10
    rmin       = 0.1
    rmax       = 10.0
    nr         = 5
    zmin       = -10.0
    zmax       = 10.0
    nz         = 10
    r0         = (rmax+rmin)/2.0
    z0         = (zmin+zmax)/2.0
    psi0       = 0
    psi1       = 1
    nrzs       = ntheta

    rs = np.cos(np.linspace(0, 2*np.math.pi, nrzs))
    zs = np.sin(np.linspace(0, 2*np.math.pi, nrzs))

    psi_rz             = np.ones((nr,nz))
    theta_psithetageom = np.ones((npsi,ntheta))
    nu_psitheta        = np.ones((npsi,ntheta))

    write_hdf5(fn, psimin, psimax, npsi, ntheta, rmin,
               rmax, nr, zmin, zmax, nz, r0, z0, psi0, psi1, psi_rz,
               theta_psithetageom, nu_psitheta, nrzs, rs, zs, desc=desc)


class Boozer(AscotData):
    """
    Object representing boozer data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
