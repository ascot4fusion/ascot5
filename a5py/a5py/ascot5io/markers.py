"""
Marker IO.
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5group import creategroup

def write_hdf5_particles(fn, N, ids, mass, charge,
                         r, phi, z, vR, vphi, vz,
                         weight, time):
    """
    Write particle marker input in hdf5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    N : int
        Number of markers
    ids : int N x 1 numpy array
        unique identifier for each marker (positive integer)
    charge : int
        charge (e)
    mass : real
        mass (amu)
    r : real N x 1 numpy array
        particle R coordinate
    phi : real N x 1 numpy array
        particle phi coordinate [deg]
    z : real N x 1 numpy array
        particle z coordinate
    vR : real N x 1 numpy array
        particle velocity R-component
    vphi : real N x 1 numpy array
        particle velocity phi-component
    vz : real N x 1 numpy array
        particle velocity z-component
    weight : real N x 1 numpy array
        particle weight (markers/s)
    time : real N x 1 numpy array
        particle initial time

    """

    mastergroup = "marker"
    subgroup    = "particle"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = creategroup(f, mastergroup, subgroup)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/n",      data=np.array([N]), dtype='i8').attrs['unit'] = '1';
        f.create_dataset(path + "/r",      data=r, dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/phi",    data=phi, dtype='f8').attrs['unit'] = 'deg';
        f.create_dataset(path + "/z",      data=z, dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/v_r",    data=vR, dtype='f8').attrs['unit'] = 'm/s';
        f.create_dataset(path + "/v_phi",  data=vphi, dtype='f8').attrs['unit'] = 'm/s';
        f.create_dataset(path + "/v_z",    data=vz, dtype='f8').attrs['unit'] = 'm/s';
        f.create_dataset(path + "/mass",   data=mass, dtype='f8').attrs['unit'] = 'amu';
        f.create_dataset(path + "/charge", data=charge, dtype='i4').attrs['unit'] = 'e';
        f.create_dataset(path + "/weight", data=weight, dtype='f8').attrs['unit'] = 'markers/s';
        f.create_dataset(path + "/time",   data=time, dtype='f8').attrs['unit'] = 's';
        f.create_dataset(path + "/id",     data=ids, dtype='i8').attrs['unit'] = '1';


def write_hdf5_guidingcenters(fn, N, ids, mass, charge,
                              r, phi, z, energy, pitch, theta,
                              weight, time):
    """
    Write guiding center marker input in hdf5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    N : int
        Number of markers
    ids : int N x 1 numpy array
        unique identifier for each marker (positive integer)
    charge : int
        charge (e)
    mass : real
        mass (amu)
    r : real N x 1 numpy array
        guiding center R coordinate
    phi : real N x 1 numpy array
        guiding center phi coordinate [deg]
    z : real N x 1 numpy array
        guiding center z coordinate
    energy : real N x 1 numpy array
        guiding center energy (eV)
    pitch : real N x 1 numpy array
        guiding center pitch (v_para/v_tot)
    theta : real N x 1 numpy array
        guiding center gyroangle (rad)
    weight : real N x 1 numpy array
        guiding center weight (markers/s)
    time : real N x 1 numpy array
        guiding center initial time

    """
    mastergroup = "marker"
    subgroup    = "guiding_center"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = creategroup(f, mastergroup, subgroup)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/n",      data=np.array([N]), dtype='i8').attrs['unit'] = '1';
        f.create_dataset(path + "/r",      data=r, dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/phi",    data=phi, dtype='f8').attrs['unit'] = 'deg';
        f.create_dataset(path + "/z",      data=z, dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/energy", data=energy, dtype='f8').attrs['unit'] = 'ev';
        f.create_dataset(path + "/pitch",  data=pitch, dtype='f8').attrs['unit'] = '1';
        f.create_dataset(path + "/theta",  data=theta, dtype='f8').attrs['unit'] = 'rad';
        f.create_dataset(path + "/mass",   data=mass, dtype='f8').attrs['unit'] = 'amu';
        f.create_dataset(path + "/charge", data=charge, dtype='i4').attrs['unit'] = 'e';
        f.create_dataset(path + "/weight", data=weight, dtype='f8').attrs['unit'] = 'markers/s';
        f.create_dataset(path + "/time",   data=time, dtype='f8').attrs['unit'] = 's';
        f.create_dataset(path + "/id",     data=ids, dtype='i8').attrs['unit'] = '1';


def write_hdf5_fieldlines(fn, N, ids, r, phi, z, pitch, weight, time):
    """
    Write magnetic field line marker input in hdf5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    N : int
        Number of markers
    ids : int N x 1 numpy array
        unique identifier for each marker (positive integer)
    r : real N x 1 numpy array
        magnetic field line R coordinate
    phi : real N x 1 numpy array
        magnetic field line phi coordinate [deg]
    z : real N x 1 numpy array
        magnetic field line z coordinate
    pitch : real N x 1 numpy array
        magnetic field line pitch whose sign defines
        the direction field line is traced
    weight : real N x 1 numpy array
        magnetic field line weight (markers/s)
    time : real N x 1 numpy array
        magnetic field line initial time

    """
    mastergroup = "marker"
    subgroup    = "field_line"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = creategroup(f, mastergroup, subgroup)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/n",      data=np.array([N]), dtype='i8').attrs['unit'] = '1';
        f.create_dataset(path + "/r",      data=r, dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/phi",    data=phi, dtype='f8').attrs['unit'] = 'deg';
        f.create_dataset(path + "/z",      data=z, dtype='f8').attrs['unit'] = 'm';
        f.create_dataset(path + "/pitch",  data=pitch, dtype='f8').attrs['unit'] = '1';
        f.create_dataset(path + "/weight", data=weight, dtype='f8').attrs['unit'] = 'markers/s';
        f.create_dataset(path + "/time",   data=time, dtype='f8').attrs['unit'] = 's';
        f.create_dataset(path + "/id",     data=ids, dtype='i8').attrs['unit'] = '1';


def read_hdf5_particles(fn, qid):
    """
    Read particle input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the particle data to be read.

    Returns
    -------

    Dictionary containing particle data.
    """

    out = {}
    with h5py.File(fn, "r") as f:
        path = "marker/particle-"+qid

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for field in f[path]:
            out[field] = f[path][field][:]

    return out

def read_hdf5_guidingcenters(fn, qid):
    """
    Read guiding-center input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the guiding-center data to be read.

    Returns
    -------

    Dictionary containing guiding-center data.
    """

    out = {}
    with h5py.File(fn, "r") as f:
        path = "marker/guiding_center-"+qid

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for field in f[path]:
            out[field] = f[path][field][:]

    return out

def read_hdf5_fieldlines(fn, qid):
    """
    Read field-line input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the field-line data to be read.

    Returns
    -------

    Dictionary containing field-line data.
    """

    out = {}
    with h5py.File(fn, "r") as f:
        path = "marker/field_line-"+qid

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for field in f[path]:
            out[field] = f[path][field][:]

    return out


def write_hdf5(fn, markers):
    # TODO move to ascot5.py
    with h5py.File(fn, "a") as f:
        if "markers" in f:
            del f["markers"]

    prt = markers["particle"]
    if prt["N"] > 0:
        write_hdf5_particles(fn, prt["N"], prt["id"], prt["mass"], prt["charge"],
                         prt["r"], prt["phi"], prt["z"], prt["vR"], prt["vphi"], prt["vz"],
                         prt["weight"], prt["time"])

    prt = markers["guiding_center"]
    if prt["N"] > 0:
        write_hdf5_guidingcenters(fn, prt["N"], prt["id"], prt["mass"], prt["charge"],
                                  prt["r"], prt["phi"], prt["z"], prt["energy"],
                                  prt["pitch"], prt["theta"],
                                  prt["weight"], prt["time"])

    prt = markers["field_line"]
    if prt["N"] > 0:
         write_hdf5_fieldlines(fn, prt["N"], prt["id"], prt["r"], prt["phi"],
                               prt["z"], prt["pitch"], prt["weight"], prt["time"])
