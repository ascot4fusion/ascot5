"""
Marker IO.

File: mrk_fl.py
"""
import h5py
import numpy as np

from . ascot5file import add_group, read_data
from a5py.ascot5io.mrk import mrk
import a5py.ascot5io.mrk

def write_hdf5(fn, n, ids, r, phi, z, pitch, weight, time, desc=None):
    """
    Write magnetic field line marker input in hdf5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        n : int <br>
            Number of markers.
        ids : array_like (n,1) <br>
            Unique identifier for each marker (must be a positive integer).
        r : array_like (n,1) <br>
            Magnetic field line R coordinate [m].
        phi : array_like (n,1) <br>
            Magnetic field line phi coordinate [deg].
        z : array_like (n,1) <br>
            Magnetic field line z coordinate [m].
        pitch : array_like (n,1) <br>
            Sign which defines the direction field line is traced, + is parallel
        weight : array_like (n,1) <br>
            Magnetic field line weight [markers/s].
        time : array_like (n,1) <br>
            Magnetic field line initial time [s].
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """
    assert ids.size    == n
    assert r.size      == n
    assert phi.size    == n
    assert z.size      == n
    assert pitch.size  == n
    assert weight.size == n
    assert time.size   == n

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
        return read_hdf5(self._file, self.get_qid())


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
