"""
"""

import copy
import numpy as np
import h5py

from ._iohelpers.treedata import DataGroup
from ._iohelpers.fileapi import read_data
from a5py.marker.plot import plot_histogram

from a5py.plotting import openfigureifnoaxes
from a5py.physlib.species import species as getspecies

class mrk(DataGroup):
    """A class acting as a superclass for all marker types.
    """

    def prune(self,keepNmarkers, data=None, probabilities=None):
        """Keep a subset of markers.

        Args:
            keepNmarkers : int <br>
                How many markers to keep.

            data=None         : dict <br>
                The data from .read() -method.
                If None, will be automatically read.

            probabilities=None
                What is the probability to include each marker
                If None, use the weights, (sum normalized to 1.0).

        Returns:
            A dictionary as if from .read()
        """

        if data is None:
            data = self.read()
        else:
            data = copy.deepcopy(data)

        n = data['n'][0]
        totalWeight = np.sum(data["weight"])

        if probabilities is None:
            probabilities = data["weight"] / totalWeight

        fortune = np.random.choice(np.arange(n),size=(keepNmarkers,),replace=False,p=probabilities)

        for k in data.keys():
            if k=='n':
                data[k][0]=keepNmarkers
                continue
            data[k] = data[k][fortune]

        newWeight = np.sum(data["weight"])

        data['weight'] *= totalWeight / newWeight

        return data

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

    def eval_rho(self, ascotpy):
        """
        Evaluate rho coordinate.
        """
        with self as h5:
            r   = read_data(h5, "r")
            phi = read_data(h5, "phi")
            z   = read_data(h5, "z")
        rho = ascotpy.evaluate(r, phi, z, 0, "rho")
        return rho

    def eval_phi(self, ascotpy):
        """
        Evaluate toroidal angle.
        """
        with self as h5:
            phi = read_data(h5, "r")

        return phi

    def eval_energy(self, ascotpy):
        """
        Evaluate energy.

        Implement in subclass.
        """
        return None

    def eval_pitch(self, ascotpy):
        """
        Evaluate pitch.

        Implemented in subclass.
        """
        return None

    @staticmethod
    def read_hdf5(fn, qid, prefix):
        """
        Read particle marker input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "marker/"+prefix+"_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                # Make all read-in datasets 1-d arrays,
                # regardless of their original dimensionality
                d=f[path][key][:]
                out[key] = np.reshape(d,newshape=(d.size,))

        out["ids"] = out["id"]
        del out["id"]
        return out

    @staticmethod
    def generatemrk(nmrk, mrktype, species=None):
        """Generate dummy marker input of given type and species.
        """
        mrk = {
            "n"      : nmrk,
            "ids"    : ( 1 + np.arange(nmrk) ),
            "r"      : np.zeros((nmrk,)),
            "z"      : np.zeros((nmrk,)),
            "phi"    : np.zeros((nmrk,)),
            "weight" : np.ones((nmrk,)),
            "time"   : np.zeros((nmrk,)),
        }
        if species is not None:
            species = getspecies(species)

        if mrktype == "particle":
            mrk["vphi"]   = np.zeros((nmrk,))
            mrk["vz"]     = np.zeros((nmrk,))
            mrk["vphi"]   = np.zeros((nmrk,))
            mrk["mass"]   = ( species["mass"]   * np.ones((nmrk,)) ).to_value("amu")
            mrk["charge"] = species["charge"] * np.ones((nmrk,), dtype=np.int16)
            mrk["anum"]   = species["anum"]   * np.ones((nmrk,), dtype=np.int16)
            mrk["znum"]   = species["znum"]   * np.ones((nmrk,), dtype=np.int16)
        if mrktype == "gc":
            mrk["pitch"]  = np.zeros((nmrk,))
            mrk["energy"] = np.zeros((nmrk,))
            mrk["zeta"]   = np.zeros((nmrk,))
            mrk["mass"]   = ( species["mass"]   * np.ones((nmrk,)) ).to_value("amu")
            mrk["charge"] = species["charge"] * np.ones((nmrk,), dtype=np.int16)
            mrk["anum"]   = species["anum"]   * np.ones((nmrk,), dtype=np.int16)
            mrk["znum"]   = species["znum"]   * np.ones((nmrk,), dtype=np.int16)
        if mrktype == "ml":
            mrk["pitch"]  = np.zeros((nmrk,))

        return mrk

class mrk_fl(mrk):
    """Object representing field line marker data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None, desc=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data, desc=desc)


    def eval_pitch(self, ascotpy):
        """
        Evaluate pitch.
        """
        with self as h5:
            pitch = read_data(h5, "pitch")

        return pitch

    @staticmethod
    def write_hdf5(fn, n, ids, r, phi, z, pitch, weight, time, desc=None):
        """Write magnetic field line marker input in hdf5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        n : int
            Number of markers.
        ids : array_like (n,1)
            Unique identifier for each marker (must be a positive integer).
        r : array_like (n,1)
            Magnetic field line R coordinate [m].
        phi : array_like (n,1)
            Magnetic field line phi coordinate [deg].
        z : array_like (n,1)
            Magnetic field line z coordinate [m].
        pitch : array_like (n,1)
            Sign which defines the direction field line is traced, + is parallel
        weight : array_like (n,1)
            Magnetic field line weight [markers/s].
        time : array_like (n,1)
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

            g.create_dataset("n",      (1,1), data=n,      dtype='i8').attrs['unit'] = '1';
            g.create_dataset("r",      (n,1), data=r,      dtype='f8').attrs['unit'] = 'm';
            g.create_dataset("phi",    (n,1), data=phi,    dtype='f8').attrs['unit'] = 'deg';
            g.create_dataset("z",      (n,1), data=z,      dtype='f8').attrs['unit'] = 'm';
            g.create_dataset("pitch",  (n,1), data=pitch,  dtype='f8').attrs['unit'] = '1';
            g.create_dataset("weight", (n,1), data=weight, dtype='f8').attrs['unit'] = 'markers/s';
            g.create_dataset("time",   (n,1), data=time,   dtype='f8').attrs['unit'] = 's';
            g.create_dataset("id",     (n,1), data=ids,    dtype='i8').attrs['unit'] = '1';

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read field line marker input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        prefix='fl'
        return a5py.ascot5io.mrk.read_hdf5(fn, qid, prefix)

class mrk_gc(mrk):
    """Object representing guiding center marker data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None, desc=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data, desc=desc)


    def eval_energy(self, ascotpy):
        with self as h5:
            energy = read_data(h5, "energy")
        return energy


    def eval_pitch(self, ascotpy):
        with self as h5:
            pitch = read_data(h5, "pitch")
        return pitch

    @staticmethod
    def write_hdf5(fn, n, ids, mass, charge, r, phi, z, energy, pitch, zeta,
                   anum, znum, weight, time, desc=None):
        """Write guiding center marker input in hdf5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        n : int
            Number of markers.
        ids : array_like (n,1)
            Unique identifier for each marker (must be a positive integer).
        charge : array_like (n,1)
            Charge [e].
        mass : array_like (n,1)
            Mass [amu].
        r : array_like (n,1)
            Guiding center R coordinate [m].
        phi : array_like (n,1)
            Guiding center phi coordinate [deg].
        z : array_like (n,1)
            Guiding center z coordinate [m].
        energy : array_like (n,1)
            Guiding center energy [eV].
        pitch : array_like (n,1)
            Guiding center pitch (v_para/v_tot).
        zeta : array_like (n,1)
            Guiding center gyroangle [rad].
        anum : array_like (n,1)
            Marker species atomic mass number.
        znum : array_like (n,1)
            Marker species charge number.
        weight : array_like (n,1)
            Guiding center weight [markers/s].
        time : array_like (n,1)
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

            g.create_dataset("n",      (1,1), data=n,      dtype='i8').attrs['unit'] = '1';
            g.create_dataset("r",      (n,1), data=r,      dtype='f8').attrs['unit'] = 'm';
            g.create_dataset("phi",    (n,1), data=phi,    dtype='f8').attrs['unit'] = 'deg';
            g.create_dataset("z",      (n,1), data=z,      dtype='f8').attrs['unit'] = 'm';
            g.create_dataset("energy", (n,1), data=energy, dtype='f8').attrs['unit'] = 'ev';
            g.create_dataset("pitch",  (n,1), data=pitch,  dtype='f8').attrs['unit'] = '1';
            g.create_dataset("zeta",   (n,1), data=zeta,   dtype='f8').attrs['unit'] = 'rad';
            g.create_dataset("mass",   (n,1), data=mass,   dtype='f8').attrs['unit'] = 'amu';
            g.create_dataset("charge", (n,1), data=charge, dtype='i4').attrs['unit'] = 'e';
            g.create_dataset("anum",   (n,1), data=anum,   dtype='i4').attrs['unit'] = '1';
            g.create_dataset("znum",   (n,1), data=znum,   dtype='i4').attrs['unit'] = '1';
            g.create_dataset("weight", (n,1), data=weight, dtype='f8').attrs['unit'] = 'markers/s';
            g.create_dataset("time",   (n,1), data=time,   dtype='f8').attrs['unit'] = 's';
            g.create_dataset("id",     (n,1), data=ids,    dtype='i8').attrs['unit'] = '1';

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read guiding center marker input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """
        prefix='gc'
        return a5py.ascot5io.mrk.read_hdf5(fn, qid, prefix)

class mrk_prt(mrk):
    """Object representing particle marker data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None, desc=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data, desc=desc)

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
    def write_hdf5(fn, n, ids, mass, charge, r, phi, z, vr, vphi, vz,
                   anum, znum, weight, time, desc=None):
        """Write particle marker input in hdf5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        n : int
            Number of markers.
        ids : array_like (n,1)
            Unique identifier for each marker (must be a positive integer).
        mass : array_like (n,1)
            Mass [amu].
        charge : array_like (n,1)
            Charge [e].
        r : array_like (n,1)
            Particle R coordinate [m].
        phi : array_like (n,1)
            Particle phi coordinate [deg].
        z : array_like (n,1)
            Particle z coordinate [m].
        vr : array_like (n,1)
            Particle velocity R-component [m/s].
        vphi : array_like (n,1)
            Particle velocity phi-component [m/s].
        vz : array_like (n,1)
            Particle velocity z-component [m/s].
        anum : array_like (n,1)
            Marker species atomic mass number.
        znum : array_like (n,1)
            Marker species charge number.
        weight : array_like (n,1)
            Particle weight [markers/s].
        time : array_like (n,1)
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

            g.create_dataset("n",      (1,1), data=n,      dtype='i8').attrs['unit'] = '1';
            g.create_dataset("r",      (n,1), data=r,      dtype='f8').attrs['unit'] = 'm';
            g.create_dataset("phi",    (n,1), data=phi,    dtype='f8').attrs['unit'] = 'deg';
            g.create_dataset("z",      (n,1), data=z,      dtype='f8').attrs['unit'] = 'm';
            g.create_dataset("vr",     (n,1), data=vr,     dtype='f8').attrs['unit'] = 'm/s';
            g.create_dataset("vphi",   (n,1), data=vphi,   dtype='f8').attrs['unit'] = 'm/s';
            g.create_dataset("vz",     (n,1), data=vz,     dtype='f8').attrs['unit'] = 'm/s';
            g.create_dataset("mass",   (n,1), data=mass,   dtype='f8').attrs['unit'] = 'amu';
            g.create_dataset("charge", (n,1), data=charge, dtype='i4').attrs['unit'] = 'e';
            g.create_dataset("anum",   (n,1), data=anum,   dtype='i4').attrs['unit'] = '1';
            g.create_dataset("znum",   (n,1), data=znum,   dtype='i4').attrs['unit'] = '1';
            g.create_dataset("weight", (n,1), data=weight, dtype='f8').attrs['unit'] = 'markers/s';
            g.create_dataset("time",   (n,1), data=time,   dtype='f8').attrs['unit'] = 's';
            g.create_dataset("id",     (n,1), data=ids,    dtype='i8').attrs['unit'] = '1';

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read particle marker input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """
        prefix='prt'
        return a5py.ascot5io.mrk.read_hdf5(fn, qid, prefix)
