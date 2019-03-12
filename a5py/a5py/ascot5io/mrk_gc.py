"""
Marker IO.

File: mrk_gc.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, n, ids, mass, charge,
               r, phi, z, energy, pitch, zeta,
               anum, znum, weight, time, desc=None):
    """
    Write guiding center marker input in hdf5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        n : int <br>
            Number of markers.
        ids : array_like (n,1) <br>
            Unique identifier for each marker (must be a positive integer).
        charge : array_like (n,1) <br>
            Charge [e].
        mass : array_like (n,1) <br>
            Mass [amu].
        r : array_like (n,1) <br>
            Guiding center R coordinate [m].
        phi : array_like (n,1) <br>
            Guiding center phi coordinate [deg].
        z : array_like (n,1) <br>
            Guiding center z coordinate [m].
        energy : array_like (n,1) <br>
            Guiding center energy [eV].
        pitch : array_like (n,1) <br>
            Guiding center pitch (v_para/v_tot).
        zeta : array_like (n,1) <br>
            Guiding center gyroangle [rad].
        anum : array_like (n,1) <br>
            Marker species atomic mass number.
        znum : array_like (n,1) <br>
            Marker species charge number.
        weight : array_like (n,1) <br>
            Guiding center weight [markers/s].
        time : array_like (n,1) <br>
            Guiding center initial time [s].
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
    assert energy.size == n
    assert pitch.size  == n
    assert zeta.size   == n
    assert anum.size   == n
    assert znum.size   == n
    assert weight.size == n
    assert time.size   == n

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

    path = "marker/mrk_gc_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


class mrk_gc(AscotData):
    """
    Object representing guiding center marker data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
