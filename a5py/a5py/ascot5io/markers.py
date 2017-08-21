"""
Marker IO.
"""
import h5py
import numpy as np
import random
import datetime

def write_hdf5_particles(fn, ids, mass, charge, 
                         r, phi, z, vR, vphi, vz, 
                         weight, time):
    """
    Write particle marker input in hdf5 file.

    TODO Not compatible with new HDF5 format.

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
    group = "markers"
    type_ = "particle"
    path = "markers/particle"

    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    if not group in f:
        o = f.create_group(group)

    # Remove group if one is already present.
    if path in f:
        del f[path]
    f.create_group(path)

    f[group].attrs['n_particle'] = np.int64_(N)

    if not "n_guiding_center" in f[group].attrs:
        f[group].attrs["n_guiding_center"] = 0

    if not "n_field_line" in f[group].attrs:
        f[group].attrs["n_field_line"] = 0

    # TODO Check that inputs are consistent.

    # Metadata.
    qid = random.getrandbits(64)
    f[path].attrs["qid"]  = np.int64_(qid)
    f[path].attrs["date"] = np.string_(datetime.datetime.now())

    # Actual data.
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

    f.close()

def write_hdf5_guidingcenters(fn, ids, mass, charge, 
                              r, phi, z, energy, pitch, theta, 
                              weight, time):
    """
    Write guiding center marker input in hdf5 file.

    TODO Not compatible with new HDF5 format.

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
    group = "markers"
    type_ = "guiding_center"
    path = "markers/guiding_center"

    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    if not group in f:
        o = f.create_group(group)

    # Remove group if one is already present.
    if path in f:
        del f[path]
    f.create_group(path)

    f[group].attrs['n_guiding_center'] = np.int64_(N)

    if not "n_particle" in f[group].attrs:
        f[group].attrs["n_particle"] = 0

    if not "n_field_line" in f[group].attrs:
        f[group].attrs["n_field_line"] = 0

    # TODO Check that inputs are consistent.

    # Metadata.
    qid = random.getrandbits(64)
    f[path].attrs["qid"]  = np.int64_(qid)
    f[path].attrs["date"] = np.string_(datetime.datetime.now())

    # Actual data.
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

    f.close()

def write_hdf5_fieldlines(fn, N, ids, r, phi, z, pitch, weight, time):
    """
    Write magnetic field line marker input in hdf5 file.

    TODO Not compatible with new HDF5 format.

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
    group = "markers"
    type_ = "field_line"
    path = "markers/field_line"

    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    if not group in f:
        o = f.create_group(group)

    # Remove group if one is already present.
    if path in f:
        del f[path]
    f.create_group(path)

    f[group].attrs['n_field_line'] = np.int64_(N)

    if not "n_particle" in f[group].attrs:
        f[group].attrs["n_particle"] = 0

    if not "n_guiding_center" in f[group].attrs:
        f[group].attrs["n_guiding_center"] = 0

    # TODO Check that inputs are consistent.

    # Metadata.
    qid = random.getrandbits(64)
    f[path].attrs["qid"]  = np.int64_(qid)
    f[path].attrs["date"] = np.string_(datetime.datetime.now())

    # Actual data.
    f.create_dataset(path + "/r",      data=r, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset(path + "/phi",    data=phi, dtype='f8').attrs['unit'] = 'deg';
    f.create_dataset(path + "/z",      data=z, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset(path + "/pitch",  data=pitch, dtype='f8').attrs['unit'] = '1';
    f.create_dataset(path + "/weight", data=weight, dtype='f8').attrs['unit'] = 'markers/s';
    f.create_dataset(path + "/time",   data=time, dtype='f8').attrs['unit'] = 's';
    f.create_dataset(path + "/id",     data=ids, dtype='i8').attrs['unit'] = '1';

    f.close()


def red_hdf5(fn):
    """
    Read marker input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing marker data.
    """
    out = {}
    f = h5py.File(fn, "r")

    out["particle"] = {}
    out["particle"]["N"] = f["markers"].attrs["n_particle"]
    out["guiding_center"] = {}
    out["guiding_center"]["N"] = f["markers"].attrs["n_guiding_center"]
    out["field_line"] = {}
    out["field_line"]["N"] = f["markers"].attrs["n_field_line"]
    
    if out["particle"]["N"] > 0:
        # Metadata.
        out["particle"]["qid"]  = f[path].attrs["qid"]
        out["particle"]["date"] = f[path].attrs["date"]

        # Actual data.
        for field in f["markers/particle"]:
            out["particle"][field] = f["markers/particle" + field][:]

    if out["guiding_center"]["N"] > 0:
        # Metadata.
        out["guiding_center"]["qid"]  = f[path].attrs["qid"]
        out["guiding_center"]["date"] = f[path].attrs["date"]

        # Actual data.
        for field in f["markers/guiding_center"]:
            out["guiding_center"][field] = f["markers/guiding_center" + field][:]

    if out["field_line"]["N"] > 0:
        # Metadata.
        out["field_line"]["qid"]  = f[path].attrs["qid"]
        out["field_line"]["date"] = f[path].attrs["date"]

        # Actual data.
        for field in f["markers/field_line"]:
            out["field_line"][field] = f["markers/field_line" + field][:]

    f.close()

    

