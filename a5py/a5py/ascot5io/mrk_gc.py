"""Marker IO.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from a5py.ascot5io.mrk import mrk
import a5py.ascot5io.mrk

class mrk_gc(mrk):
    """
    Object representing guiding center marker data.
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
