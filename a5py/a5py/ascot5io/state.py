"""
State HDF5 IO module.

File: state.py
"""

import numpy as np
import h5py

from a5py.marker.alias import get as alias
import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.plot as plot
import a5py.marker.endcond as endcondmod

from a5py.marker.alias import get as alias

from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
               b_phimin, b_phimax, b_nphi,
               axisr, axisz, psi, psi0, psi1, br, bphi, bz,
               psi_rmin=None, psi_rmax=None, psi_nr=None,
               psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
    dict_keys(['anum', 'anum_unit',
               'bphi', 'bphi_unit', 'bphidphi', 'bphidphi_unit',
               'bphidr', 'bphidr_unit', 'bphidz', 'bphidz_unit',
               'br', 'br_unit', 'brdphi', 'brdphi_unit',
               'brdr', 'brdr_unit', 'brdz', 'brdz_unit',
               'bz', 'bz_unit', 'bzdphi', 'bzdphi_unit',
               'bzdr', 'bzdr_unit', 'bzdz', 'bzdz_unit',
               'charge', 'charge_unit', 'cputime', 'cputime_unit',
               'endcond', 'endcond_unit', 'errorline', 'errorline_unit',
               'errormod', 'errormod_unit', 'errormsg', 'errormsg_unit',
               'id', 'id_unit', 'mass', 'mass_unit', 'mu', 'mu_unit',
               'phi', 'phi_unit', 'phiprt', 'phiprt_unit',
               'r', 'r_unit', 'rprt', 'rprt_unit',
               'rho', 'rho_unit', 'theta', 'theta_unit',
               'time', 'time_unit', 'vpar', 'vpar_unit',
               'vphi', 'vphi_unit', 'vr', 'vr_unit',
               'vz', 'vz_unit', 'walltile', 'walltile_unit',
               'weight', 'weight_unit',
               'z', 'z_unit', 'zeta', 'zeta_unit',
               'znum', 'znum_unit', 'zprt', 'zprt_unit', 'N'])
    """
    Write 3DS magnetic field input in HDF5 file.

    Note that br and bz should not include the equilibrium component of the
    magnetic field as that is calculated from psi by ASCOT5 during the
    simulation.

    It is possible to use different Rz grids for psi and magnetic field
    components by giving Rz grid for psi separately.

    The toroidal angle phi is treated as a periodic coordinate meaning that
    B(phi) = B(phi + n*(b_phimax - b_phimin)). Do note that to avoid dublicate
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
        bphi : array_like (nr,nphi,nz) <br>
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nr,nphi,nz) <br>
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
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
    """

    parent = "bfield"
    group  = "B_3DS"
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

    assert psi.shape  == (psi_nr,psi_nz)
    assert br.shape   == (b_nr,b_nphi,b_nz)
    assert bphi.shape == (b_nr,b_nphi,b_nz)
    assert bz.shape   == (b_nr,b_nphi,b_nz)

    psi  = np.transpose(psi)
    br   = np.transpose(br,   (2,1,0))
    bphi = np.transpose(bphi, (2,1,0))
    bz   = np.transpose(bz,   (2,1,0))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("b_rmin",   (1,), data=b_rmin,   dtype="f8")
        g.create_dataset("b_rmax",   (1,), data=b_rmax,   dtype="f8")
        g.create_dataset("b_nr",     (1,), data=b_nr,     dtype="i4")
        g.create_dataset("b_phimin", (1,), data=b_phimin, dtype="f8")
        g.create_dataset("b_phimax", (1,), data=b_phimax, dtype="f8")
        g.create_dataset("b_nphi",   (1,), data=b_nphi,   dtype="i4")
        g.create_dataset("b_zmin",   (1,), data=b_zmin,   dtype="f8")
        g.create_dataset("b_zmax",   (1,), data=b_zmax,   dtype="f8")
        g.create_dataset("b_nz",     (1,), data=b_nz,     dtype="i4")
        g.create_dataset("psi_rmin", (1,), data=psi_rmin, dtype="f8")
        g.create_dataset("psi_rmax", (1,), data=psi_rmax, dtype="f8")
        g.create_dataset("psi_nr",   (1,), data=psi_nr,   dtype="i4")
        g.create_dataset("psi_zmin", (1,), data=psi_zmin, dtype="f8")
        g.create_dataset("psi_zmax", (1,), data=psi_zmax, dtype="f8")
        g.create_dataset("psi_nz",   (1,), data=psi_nz,   dtype="i4")
        g.create_dataset("axisr",    (1,), data=axisr,    dtype="f8")
        g.create_dataset("axisz",    (1,), data=axisz,    dtype="f8")
        g.create_dataset("psi0",     (1,), data=psi0,     dtype="f8")
        g.create_dataset("psi1",     (1,), data=psi1,     dtype="f8")

        g.create_dataset("psi",  (psi_nz, psi_nr),     data=psi,  dtype="f8")
        g.create_dataset("br",   (b_nz, b_nphi, b_nr), data=br,   dtype="f8")
        g.create_dataset("bphi", (b_nz, b_nphi, b_nr), data=bphi, dtype="f8")
        g.create_dataset("bz",   (b_nz, b_nphi, b_nr), data=bz,   dtype="f8")

    return gname

def read_hdf5(fn, qid, name):
    """
    Read all or specified states.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    read : string list, optional
        Which states are read e.g. "inistate". Default is all.

    Returns
    -------

    Dictionary storing the states that were read.
    """

    path = "results/run_" + qid + "/" + name

    with h5py.File(fn,"r") as f:
        out = {}

        for field in f[path]:
            out[field]           = f[path][field][:]
            out[field + "_unit"] = f[path][field].attrs["unit"]

        # TODO Parse endconditions.

        # Find number of markers and check that no markers share same id
        # (which they shouldn't).
        out["N"] = np.unique(out["id"]).size
        if out["N"] != out["id"].size:
            print("Warning: Markers don't have unique Id.")

    return out


class State(AscotData):
    """
    """

    def __init__(self, hdf5, runnode):
        """
        Initialize state object from given HDF5 file to given RunNode.
        """
        self._runnode = runnode
        super().__init__(hdf5)


    def read(self):
        """
        Read orbit data to dictionary.
        """
        return read_hdf5(self._file, self.get_qid(), self._path.split("/")[-1])


    def __getitem__(self, key):
        """
        Return queried quantity.

        The quantity is returned as a single numpy array ordered by id and time.
        Internally, this method first sees if the quantity can be read directly
        from HDF5. If not, then it tries to see if it is present in endstate and
        can be copied from there (e.g. mass). If not, then quantity is evaluated
        by first determining if the stored orbit type is field line (has no
        charge), guiding center (has magnetic moment), or particle.

        Args:
            key : str <br>
                Name of the quantity (see alias.py for a complete list).
        Returns:
            The quantity in SI units ordered by id and time.
        """
        with self as h5:

            key  = alias(key)
            item = None

            # See if the field can be read directly and without conversions
            h5keys = list(h5.keys())
            if key in h5keys:
                item = h5[key][:]

                # Unit conversions
                if key == "charge":
                    f    = lambda x: interpret.charge_C(x)
                    item = np.array([f(x) for x in item]).ravel()
                if key == "mu":
                    f    = lambda x: interpret.energy_J(x)
                    item = np.array([f(x) for x in item]).ravel()
                if key == "phi":
                    item = item * np.pi/180
                if key == "mass":
                    f    = lambda x: interpret.mass_kg(x)
                    item = np.array([f(x) for x in item]).ravel()
                if key == "endcond":
                    err = h5["errormsg"][:]
                    item = item << 2
                    item[err > 0] = item[err > 0] & endcondmod.getbin("aborted")
                    item[item==0] = endcondmod.getbin("none")

            if item is None:

                # Convert guiding-center quantities to SI units
                f      = lambda x: interpret.mass_kg(x)
                mass   = np.array([f(x) for x in h5["mass"][:]]).ravel()
                f      = lambda x: interpret.charge_C(x)
                charge = np.array([f(x) for x in h5["charge"][:]]).ravel()
                f      = lambda x: interpret.energy_J(x)
                mu     = np.array([f(x) for x in h5["mu"][:]]).ravel()
                phi    = h5["phi"][:] * np.pi/180

                try:
                    item = marker.eval_guidingcenter(
                        key, mass=mass, charge=charge,
                        R=h5["r"][:], phi=phi, z=h5["z"][:],
                        mu=mu, vpar=h5["vpar"][:],
                        theta=h5["theta"][:],
                        BR=h5["br"][:], Bphi=h5["bphi"][:],
                        Bz=h5["bz"][:])
                except ValueError:
                    pass

                if item is None:
                    item = marker.eval_particle(
                        key, mass=mass, charge=charge,
                        R=h5["r"][:], phi=phi, z=h5["z"][:],
                        vR=h5["vr"][:], vphi=h5["vphi"][:], vz=h5["vz"][:],
                        BR=h5["br"][:], Bphi=h5["bphi"][:], Bz=h5["bz"][:])

            # Order by id
            ids  = h5["id"][:]
            time = h5["time"][:]
            idx  = np.lexsort((time, ids))

            return item[idx]


    def get(self, key, ids=None, endcond=None, pncrid=None, SI=True):
        """
        Same as __getitem__ but with option to filter which points are returned.

        Args:
            key : str <br>
                Name of the quantity (see alias.py for a complete list).
            ids : int, array_like, optional <br>
                Id or a list of ids whose data points are returned.
            endcond : str, array_like, optional <br>
                Endcond of those  markers which are returned.
            pncrid : str, array_like, optional <br>
                Poincare ID of those  markers which are returned.
            SI : bool, optional <br>
                Whether to return quantity in SI units or Ascot units.
        Returns:
            The quantity.
        """
        val = self[key]

        idx = np.ones(val.shape, dtype=bool)

        if endcond is not None:
            if hasattr(self._runnode, "endstate"):
                ec = self._runnode.endstate["endcond"]
            else:
                ec = self["endcond"]

            idx = np.logical_and( idx, ec == endcondmod.getbin(endcond) )

        if pncrid is not None:
            idx = np.logical_and(idx, self["pncrid"] == pncrid)

        val = val[idx]

        if not SI:
            key = alias(key)

            if key in ["energy", "mu"]:
                f      = lambda x: interpret.energy_eV(x)
                val   = np.array([f(x) for x in val]).ravel()

            if key in ["phi", "phimod"]:
                val = val*180/np.pi

            if key in ["mass"]:
                f      = lambda x: interpret.mass_amu(x)
                val   = np.array([f(x) for x in val]).ravel()

            if key in ["charge"]:
                f      = lambda x: interpret.charge_e(x)
                val   = np.array([f(x) for x in val]).ravel()

        return val


    def listendconds(self):
        ec, counts = np.unique(self["endcond"], return_counts=True)
        endcond = []
        for e in ec:
            endcond.append(endcondmod.getname(e))

        return (endcond, counts)


    def scatter(self, x=None, y=None, z=None, c=None, endcond=None, pncrid=None,
                equal=False, log=False, axes=None):
        """
        Make scatter plot.

        """
        ids = self.get("id", endcond=endcond, pncrid=pncrid)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, endcond=endcond, pncrid=pncrid, SI=False)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, pncrid=pncrid, SI=False)

        zc = None
        if z is not None:
            zc = self.get(z, endcond=endcond, pncrid=pncrid, SI=False)

        cc = None
        if c is not None:
            cc = self.get(c, endcond=endcond, pncrid=pncrid, SI=False)

        if isinstance(log, tuple):
            if log[0]:
                xc = np.log10(np.absolute(xc))
            if log[1]:
                yc = np.log10(np.absolute(yc))
            if z is not None and log[2]:
                zc = np.log10(np.absolute(zc))
            if c is not None and log[2]:
                cc = np.log10(np.absolute(cc))
        elif log:
            xc = np.log10(np.absolute(xc))
            yc = np.log10(np.absolute(yc))
            if z is not None:
                zc = np.log10(np.absolute(zc))
            if c is not None:
                cc = np.log10(np.absolute(cc))

        plot.plot_scatter(x=xc, y=yc, z=zc, c=cc, equal=equal,
                          xlabel=x, ylabel=y, zlabel=z, axes=axes)


    def histogram(self, x, y=None, xbins=None, ybins=None, weight=False,
                  logx=False, logy=False, logscale=False, endcond=None,
                  axes=None):
        """
        Make histogram plot.
        """
        if y is None:
            # 1D plot

            if logx:
                xbins = np.linspace(np.log10(xbins[0]), np.log10(xbins[1]),
                                    xbins[2])
            else:
                xbins = np.linspace(xbins[0], xbins[1], xbins[2])

            weights=None
            if endcond is not None or not hasattr(self._runnode, "endstate"):
                xc = self.get(x, endcond=endcond, SI=False)
                if weight:
                    weights = self.get("weight", endcond=endcond)
                if logx:
                    xc = np.log10(np.absolute(xc))
            else:
                xc = []
                if weight:
                    weights = []

                ecs, count = self._runnode.endstate.listendconds()
                for ec in ecs:
                    xc0 = self.get(x, endcond=ec, SI=False)
                    if weight:
                        weights.append(self.get("weight", endcond=ec))
                    if logx:
                        xc0 = np.log10(np.absolute(xc0))

                    xc.append(xc0)

            plot.plot_histogram(x=xc, y=None, xbins=xbins, ybins=ybins,
                                weights=weights, logscale=logscale,
                                xlabel=x, ylabel=y, axes=axes)

        else:
            # 2D plot
            xc = self.get(x, endcond=endcond, SI=False)
            yc = self.get(y, endcond=endcond, SI=False)

            if logx:
                xc = np.log10(np.absolute(xc))
                xbins = np.linspace(np.log10(xbins[0]), np.log10(xbins[1]),
                                    xbins[2])
            else:
                xbins = np.linspace(xbins[0], xbins[1], xbins[2])

            if logy:
                yc = np.log10(np.absolute(yc))
                ybins = np.linspace(np.log10(ybins[0]), np.log10(ybins[1]),
                                    ybins[2])
            else:
                ybins = np.linspace(ybins[0], ybins[1], ybins[2])

            weights=None
            if weight:
                weights = self.get("weight", endcond=endcond)

            plot.plot_histogram(x=xc, y=yc, xbins=xbins, ybins=ybins,
                                weights=weights, logscale=logscale,
                                xlabel=x, ylabel=y, axes=axes)
