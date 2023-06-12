"""Marker IO.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group, read_data
from a5py.ascot5io.mrk import mrk
import a5py.ascot5io.mrk

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


class mrk_fl(mrk):
    """
    Object representing field line marker data.
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
