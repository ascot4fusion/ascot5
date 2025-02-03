"""Markers whose orbits are solved within the ASCOT5 code.
"""
import copy
import numpy as np
import h5py
import unyt


class DataGroup():
    pass

from a5py.routines.plotting import openfigureifnoaxes
from a5py.physlib import parseunits
from a5py.physlib.species import species as getspecies

class Marker(DataGroup):
    """A class acting as a superclass for all marker types.
    """

    def prune(self, pick, normalize=True, mrk=None):
        """Take a subset of markers.

        Parameters
        ----------
        pick : int or array_like
            Either number of markers to be picked (in which case they are chosen
            randomly) or the marker IDs.
        normalize : bool, optional
            Reweight the output markers so that the total weight is
            the same as in input.
        mrk : dict, optional
            Marker data dictionary (same format as in corresponding
            ``write_hdf5`` function).

            If None, the data is read from the file.

        Returns
        -------
        out : dict
            Marker data dictionary (same format as in corresponding
            ``write_hdf5`` function) with picked markers.
        """

        if mrk is None:
            mrk = self.read()
        else:
            mrk = copy.deepcopy(mrk)

        weighttotal = np.sum(mrk["weight"])

        if isinstance(pick, int):
            rng  = np.random.default_rng()
            pick = rng.choice(mrk["ids"].ravel(), size=(pick,), replace=False)

        idx = np.where(np.in1d(mrk["ids"], pick))
        for k in mrk:
            if k == "n": continue
            mrk[k] = mrk[k][idx]

        if normalize:
            mrk["weight"][:] *= weighttotal / np.sum(mrk["weight"])

        mrk["n"] = pick.size
        return mrk

    @openfigureifnoaxes(projection=None)
    def plot_hist_rhophi(self, ascotpy, rbins=10, pbins=10, weighted=False,
                         axes=None):
        """
        Plot marker rho-phi histogram
        """
        weights = None
        if weighted:
            weights = self.read()["weight"]

        rho = self.eval_rho(ascotpy)
        phi = self.eval_phi(ascotpy)
        plot_histogram(rho, xbins=rbins, y=phi, ybins=pbins,
                       weights=weights,
                       logscale=False, xlabel="Normalized poloidal flux",
                       ylabel="Toroidal angle [deg]",
                       axes=axes)

    @openfigureifnoaxes(projection=None)
    def plot_hist_energypitch(self, ascotpy, ebins=10, pbins=10, weighted=False,
                              axes=None):
        """
        Plot marker energy-pitch histogram

        Energy is on logscale.
        """
        weights = None
        if weighted:
            weights = self.read()["weight"]

        pitch  = self.eval_pitch(ascotpy)
        energy = self.eval_energy(ascotpy)
        if energy is None:
            # Field lines don't have energy
            plot_histogram(x=pitch, xbins=pbins, weights=weights,
                           logscale=False, xlabel=r"Pitch [$p_\parallel/p$]",
                           axes=axes)
        else:
            energy = np.log10(energy.to("eV"))
            plot_histogram(energy, xbins=ebins, y=pitch, ybins=pbins,
                           weights=weights,
                           logscale=False, xlabel="Energy [eV]",
                           ylabel=r"Pitch [$p_\parallel/p$]",
                           axes=axes)

    @staticmethod
    def generate(mrktype, n, species=None):
        """Generate dummy marker input of given type and species.

        Parameters
        ----------
        mrktype : {"prt", "gc", "fl"}
            Type of marker input to be created.
        n : int
            Number of markers.
        species : str, optional
            Marker species to initialize anum, znum, mass, and charge.

        Returns
        -------
        mrk : dict
            Markers parameters that can be supplied to :meth:`Prt.write_hdf5`,
            :meth:`GC.write_hdf5` or :meth:`FL.write_hdf5` depending on
            ``mrktype`` value.
        """
        mrk = {
            "n"      : n,
            "ids"    : ( 1 + np.arange(n) ),
            "r"      : np.zeros((n,))*unyt.m,
            "z"      : np.zeros((n,))*unyt.m,
            "phi"    : np.zeros((n,))*unyt.deg,
            "weight" : np.ones((n,))*unyt.particles/unyt.s,
            "time"   : np.zeros((n,))*unyt.s,
        }
        if species is None:
            species = {"mass":1*unyt.amu, "charge":1*unyt.e, "anum":1, "znum":1}
        else:
            species = getspecies(species)

        if mrktype == "prt":
            mrk["vr"]     = np.zeros((n,)) * unyt.m/unyt.s
            mrk["vz"]     = np.zeros((n,)) * unyt.m/unyt.s
            mrk["vphi"]   = np.zeros((n,)) * unyt.m/unyt.s
            mrk["mass"]   = species["mass"] * np.ones((n,))
            mrk["charge"] = species["charge"] * np.ones((n,), dtype=np.int16)
            mrk["anum"]   = species["anum"]   * np.ones((n,), dtype=np.int16)
            mrk["znum"]   = species["znum"]   * np.ones((n,), dtype=np.int16)
        if mrktype == "gc":
            mrk["pitch"]  = np.zeros((n,)) * unyt.dimensionless
            mrk["energy"] = np.zeros((n,)) * unyt.eV
            mrk["zeta"]   = np.zeros((n,)) * unyt.rad
            mrk["mass"]   = species["mass"] * np.ones((n,))
            mrk["charge"] = species["charge"] * np.ones((n,), dtype=np.int16)
            mrk["anum"]   = species["anum"]   * np.ones((n,), dtype=np.int16)
            mrk["znum"]   = species["znum"]   * np.ones((n,), dtype=np.int16)
        if mrktype == "fl":
            mrk["pitch"]  = np.zeros((n,)) * unyt.dimensionless

        return mrk

class FL(Marker):
    """Magnetic-field-line tracer input.

    These are markers that trace magnetic field lines exactly. They don't
    represent physical particles except that for calculating quantities like
    ``mileage``, they are assumed to be photons travelling at the speed of
    light.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]
                if key == "n":
                    out[key] = int(out[key])

        out["ids"] = out["id"]
        del out["id"]
        return out

    def eval_pitch(self, ascotpy):
        """
        Evaluate pitch.
        """
        with self as h5:
            pitch = read_data(h5, "pitch")

        return pitch

    @staticmethod
    def write_hdf5(fn, n, ids, r, phi, z, pitch, weight, time, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        n : int
            Number of markers.
        ids : array_like (n,)
            Unique identifier for each marker (must be a positive integer).
        r : array_like (n,)
            Magnetic field line R coordinate [m].
        phi : array_like (n,)
            Magnetic field line phi coordinate [deg].
        z : array_like (n,)
            Magnetic field line z coordinate [m].
        pitch : array_like (n,)
            Sign which defines the direction field line is traced, + is parallel
        weight : array_like (n,)
            Magnetic field line weight [markers/s].
        time : array_like (n,)
            Magnetic field line initial time [s].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if ids.size    != n: raise ValueError("Inconsistent size for ids.")
        if r.size      != n: raise ValueError("Inconsistent size for r.")
        if phi.size    != n: raise ValueError("Inconsistent size for phi.")
        if z.size      != n: raise ValueError("Inconsistent size for z.")
        if pitch.size  != n: raise ValueError("Inconsistent size for pitch.")
        if weight.size != n: raise ValueError("Inconsistent size for weight.")
        if time.size   != n: raise ValueError("Inconsistent size for time.")

        parent = "marker"
        group  = "fl"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            write_data(g, "n",      n,      (1,), "i8",                compress=False)
            write_data(g, "r",      r,      (n,), "f8", "m",           compress=True )
            write_data(g, "phi",    phi,    (n,), "f8", "deg",         compress=True )
            write_data(g, "z",      z,      (n,), "f8", "m",           compress=True )
            write_data(g, "pitch",  pitch,  (n,), "f8", "1",           compress=True )
            write_data(g, "weight", weight, (n,), "f8", "particles/s", compress=True )
            write_data(g, "time",   time,   (n,), "f8", "s",           compress=True )
            write_data(g, "id",     ids,    (n,), "i8", "1",           compress=True )

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return Marker.generate(n=1, mrktype="fl", species=None)

class GC(Marker):
    """Particle input in guiding-center coordinates.

    The guiding center phase-space is (R, Phi, Z, mu, ppar, zeta), where
    R, Phi, Z are the guiding center position in cylindrical coordinates,
    mu is the magnetic moment (first adiabatic invariant), and ppar is
    the parallel momentum along the magnetic field line). The gyro-angle
    zeta is an ignorable quantity in guiding-center integration, but it is
    used in guiding-center-to-particle transformation.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]
                if key == "n":
                    out[key] = int(out[key])

        out["ids"] = out["id"]
        del out["id"]
        return out

    def eval_energy(self, ascotpy):
        with self as h5:
            energy = read_data(h5, "energy")
        return energy


    def eval_pitch(self, ascotpy):
        with self as h5:
            pitch = read_data(h5, "pitch")
        return pitch

    @staticmethod
    @parseunits(strip=True, mass="amu", charge="e", r="m", phi="deg", z="m",
                energy="eV", pitch="1", zeta="rad", anum="1", znum="1",
                time="s", weight="1/s")
    def write_hdf5(fn, n, ids, mass, charge, r, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        n : int
            Number of markers.
        ids : array_like (n,)
            Unique identifier for each marker (must be a positive integer).
        mass : array_like (n,)
            Mass [amu].
        charge : array_like (n,)
            Charge [e].
        r : array_like (n,)
            Guiding center R coordinate [m].
        phi : array_like (n,)
            Guiding center phi coordinate [deg].
        z : array_like (n,)
            Guiding center z coordinate [m].
        energy : array_like (n,)
            Guiding center energy [eV].
        pitch : array_like (n,)
            Guiding center pitch (v_para/v_tot).
        zeta : array_like (n,)
            Guiding center gyroangle [rad].
        anum : array_like (n,)
            Marker species atomic mass number.
        znum : array_like (n,)
            Marker species charge number.
        weight : array_like (n,)
            Guiding center weight [markers/s].
        time : array_like (n,)
            Guiding center initial time [s].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if ids.size    != n: raise ValueError("Inconsistent size for ids.")
        if mass.size   != n: raise ValueError("Inconsistent size for mass.")
        if charge.size != n: raise ValueError("Inconsistent size for charge.")
        if r.size      != n: raise ValueError("Inconsistent size for r.")
        if phi.size    != n: raise ValueError("Inconsistent size for phi.")
        if z.size      != n: raise ValueError("Inconsistent size for z.")
        if energy.size != n: raise ValueError("Inconsistent size for energy.")
        if pitch.size  != n: raise ValueError("Inconsistent size for pitch.")
        if zeta.size   != n: raise ValueError("Inconsistent size for zeta.")
        if anum.size   != n: raise ValueError("Inconsistent size for anum.")
        if znum.size   != n: raise ValueError("Inconsistent size for znum.")
        if weight.size != n: raise ValueError("Inconsistent size for weight.")
        if time.size   != n: raise ValueError("Inconsistent size for time.")

        parent = "marker"
        group  = "gc"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("n",      (1,), data=n,      dtype='i8').attrs['unit'] = '1';
            g.create_dataset("r",      (n,), data=r,      dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm';
            g.create_dataset("phi",    (n,), data=phi,    dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'deg';
            g.create_dataset("z",      (n,), data=z,      dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm';
            g.create_dataset("energy", (n,), data=energy, dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'ev';
            g.create_dataset("pitch",  (n,), data=pitch,  dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = '1';
            g.create_dataset("zeta",   (n,), data=zeta,   dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'rad';
            g.create_dataset("mass",   (n,), data=mass,   dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'amu';
            g.create_dataset("charge", (n,), data=charge, dtype='i4', compression="gzip", compression_opts=9).attrs['unit'] = 'e';
            g.create_dataset("anum",   (n,), data=anum,   dtype='i4', compression="gzip", compression_opts=9).attrs['unit'] = '1';
            g.create_dataset("znum",   (n,), data=znum,   dtype='i4', compression="gzip", compression_opts=9).attrs['unit'] = '1';
            g.create_dataset("weight", (n,), data=weight, dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'markers/s';
            g.create_dataset("time",   (n,), data=time,   dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 's';
            g.create_dataset("id",     (n,), data=ids,    dtype='i8', compression="gzip", compression_opts=9).attrs['unit'] = '1';

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return Marker.generate(n=1, mrktype="gc", species="alpha")

class Prt(Marker):
    """Marker input representing physical particles.

    Particle phase-space is (r, phi, z, pr, pphi, pz), where r, phi, and z are
    particle position in cylindrical coordinates, and pr, phi, and pz are its
    momentum components.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]
                if key == "n":
                    out[key] = int(out[key])

        out["ids"] = out["id"]
        del out["id"]
        return out

    def eval_energy(self, ascotpy):
        with self as h5:
            vr   = read_data(h5, "vr")
            vz   = read_data(h5, "vz")
            vphi = read_data(h5, "vphi")
            mass = read_data(h5, "mass")

        v = np.sqrt(vr**2 + vz**2 + vphi**2)
        return energy_velocity(mass, v)

    def eval_pitch(self, ascotpy):
        with self as h5:
            r    = read_data(h5, "r")
            z    = read_data(h5, "z")
            phi  = read_data(h5, "phi")
            vr   = read_data(h5, "vr")
            vz   = read_data(h5, "vz")
            vphi = read_data(h5, "vphi")
            mass = read_data(h5, "mass")

        br    = ascotpy.evaluate(r, phi, z, 0, "br")
        bz    = ascotpy.evaluate(r, phi, z, 0, "bz")
        bphi  = ascotpy.evaluate(r, phi, z, 0, "bphi")
        b     = np.sqrt(br**2 + bz**2 + bphi**2)
        v     = np.sqrt(vr**2 + vz**2 + vphi**2)
        pitch = ( vr * br + vz * bz + vphi * bphi ) / ( b * v )
        return pitch

    @staticmethod
    @parseunits(strip=True, mass="amu", charge="e", r="m", phi="deg", z="m",
                vr="m/s", vphi="m/s", vz="m/s", anum="1", znum="1",
                time="s", weight="1/s")
    def write_hdf5(fn, n, ids, mass, charge, r, phi, z, vr, vphi, vz,
                   anum, znum, weight, time, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        n : int
            Number of markers.
        ids : array_like (n,)
            Unique identifier for each marker (must be a positive integer).
        mass : array_like (n,)
            Mass [amu].
        charge : array_like (n,)
            Charge [e].
        r : array_like (n,)
            Particle R coordinate [m].
        phi : array_like (n,)
            Particle phi coordinate [deg].
        z : array_like (n,)
            Particle z coordinate [m].
        vr : array_like (n,)
            Particle velocity R-component [m/s].
        vphi : array_like (n,)
            Particle velocity phi-component [m/s].
        vz : array_like (n,)
            Particle velocity z-component [m/s].
        anum : array_like (n,)
            Marker species atomic mass number.
        znum : array_like (n,)
            Marker species charge number.
        weight : array_like (n,)
            Particle weight [markers/s].
        time : array_like (n,)
            Particle initial time [s].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if ids.size    != n: raise ValueError("Inconsistent size for ids.")
        if mass.size   != n: raise ValueError("Inconsistent size for mass.")
        if charge.size != n: raise ValueError("Inconsistent size for charge.")
        if r.size      != n: raise ValueError("Inconsistent size for r.")
        if phi.size    != n: raise ValueError("Inconsistent size for phi.")
        if z.size      != n: raise ValueError("Inconsistent size for z.")
        if vr.size     != n: raise ValueError("Inconsistent size for vR.")
        if vphi.size   != n: raise ValueError("Inconsistent size for vphi.")
        if vz.size     != n: raise ValueError("Inconsistent size for vz.")
        if anum.size   != n: raise ValueError("Inconsistent size for anum.")
        if znum.size   != n: raise ValueError("Inconsistent size for znum.")
        if weight.size != n: raise ValueError("Inconsistent size for weight.")
        if time.size   != n: raise ValueError("Inconsistent size for time.")

        parent = "marker"
        group  = "prt"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("n",      (1,), data=n,      dtype='i8').attrs['unit'] = '1';
            g.create_dataset("r",      (n,), data=r,      dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm';
            g.create_dataset("phi",    (n,), data=phi,    dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'deg';
            g.create_dataset("z",      (n,), data=z,      dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm';
            g.create_dataset("vr",     (n,), data=vr,     dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm/s';
            g.create_dataset("vphi",   (n,), data=vphi,   dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm/s';
            g.create_dataset("vz",     (n,), data=vz,     dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'm/s';
            g.create_dataset("mass",   (n,), data=mass,   dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'amu';
            g.create_dataset("charge", (n,), data=charge, dtype='i4', compression="gzip", compression_opts=9).attrs['unit'] = 'e';
            g.create_dataset("anum",   (n,), data=anum,   dtype='i4', compression="gzip", compression_opts=9).attrs['unit'] = '1';
            g.create_dataset("znum",   (n,), data=znum,   dtype='i4', compression="gzip", compression_opts=9).attrs['unit'] = '1';
            g.create_dataset("weight", (n,), data=weight, dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 'markers/s';
            g.create_dataset("time",   (n,), data=time,   dtype='f8', compression="gzip", compression_opts=9).attrs['unit'] = 's';
            g.create_dataset("id",     (n,), data=ids,    dtype='i8', compression="gzip", compression_opts=9).attrs['unit'] = '1';

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return Marker.generate(n=1, mrktype="prt", species="alpha")
