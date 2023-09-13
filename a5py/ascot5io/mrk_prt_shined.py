"""
Marker IO.

File: mrk_prt_shined.py
"""
import h5py
import numpy as np

from .coreio.fileapi import add_group
from a5py.ascot5io.mrk import mrk
import a5py.ascot5io.mrk
import a5py.ascot5io.mrk_shined


def write_hdf5(fn, n, ids, mass, charge,
               r, phi, z, vr, vphi, vz,
               anum, znum, weight, time, desc=None):
    """
    Write particle marker input in hdf5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        n : int <br>
            Number of markers.
        ids : array_like (n,1) <br>
            Unique identifier for each marker (must be a positive integer).
        mass : array_like (n,1) <br>
            Mass [amu].
        charge : array_like (n,1) <br>
            Charge [e].
        r : array_like (n,1) <br>
            Particle R coordinate [m].
        phi : array_like (n,1) <br>
            Particle phi coordinate [deg].
        z : array_like (n,1) <br>
            Particle z coordinate [m].
        vr : array_like (n,1) <br>
            Particle velocity R-component [m/s].
        vphi : array_like (n,1) <br>
            Particle velocity phi-component [m/s].
        vz : array_like (n,1) <br>
            Particle velocity z-component [m/s].
        anum : array_like (n,1) <br>
            Marker species atomic mass number.
        znum : array_like (n,1) <br>
            Marker species charge number.
        weight : array_like (n,1) <br>
            Particle weight [markers/s].
        time : array_like (n,1) <br>
            Particle initial time [s].
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """
    assert ids.size    == n
    assert mass.size   == n
    assert charge.size == n
    assert r.size      == n
    assert phi.size    == n
    assert z.size      == n
    assert vr.size     == n
    assert vphi.size   == n
    assert vz.size     == n
    assert anum.size   == n
    assert znum.size   == n
    assert weight.size == n
    assert time.size   == n

    parent = "marker_shined"
    group  = "prt_shined"
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
    prefix='prt_shined'
    return a5py.ascot5io.mrk_shined.read_hdf5(fn, qid, prefix)


class mrk_prt_shined(mrk):
    """
    Object representing particle marker data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None, desc=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data, desc=desc)
