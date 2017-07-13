import h5py
import numpy as np

def write_hdf5_particles(fn, ids, mass, charge, r, phi, z, vr, vphi, vz, weight, time):
    """Write particle marker input in hdf5 file.

    Keyword arguments:
    fn     -- path to hdf5 file
    ids    -- unique identifier for each marker (positive integer)
    mass   -- mass (amu)
    charge -- charge (e)
    r      -- particle R coordinate (m)
    phi    -- particle phi coordinate (deg)
    z      -- particle z coordinate (m)
    vr     -- particle R velocity (m/s)
    vphi   -- particle phi velocity (m/s)
    vz     -- particle z velocity (m/s)
    weight -- particle weight (markers/s)
    time   -- particle initial time (s)

    All particle inputs must be in numpy array format.
    """
    f = h5py.File(fn, "a")
    if not "/markers" in f:
        f.create_group('markers')
    
    if  "/markers/particle" in f:
        del f["/markers/particle"]
        

    f.create_group('markers/particle')

    f['markers'].attrs['n_particle'] = np.size(ids)

    if not 'n_guiding_center' in f['markers'].attrs:
        f['markers'].attrs['n_guiding_center'] = 0

    if not 'n_field_line' in f['markers'].attrs:
        f['markers'].attrs['n_field_line'] = 0

    f.create_dataset('markers/particle/r', data=r, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset('markers/particle/phi', data=phi, dtype='f8').attrs['unit'] = 'deg';
    f.create_dataset('markers/particle/z', data=z, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset('markers/particle/v_r', data=vr, dtype='f8').attrs['unit'] = 'm/s';
    f.create_dataset('markers/particle/v_phi', data=vphi, dtype='f8').attrs['unit'] = 'm/s';
    f.create_dataset('markers/particle/v_z', data=vz, dtype='f8').attrs['unit'] = 'm/s';
    f.create_dataset('markers/particle/mass', data=mass, dtype='f8').attrs['unit'] = 'amu';
    f.create_dataset('markers/particle/charge', data=charge, dtype='i4').attrs['unit'] = 'e';
    f.create_dataset('markers/particle/weight', data=weight, dtype='f8').attrs['unit'] = 'markers/s';
    f.create_dataset('markers/particle/time', data=time, dtype='f8').attrs['unit'] = 's';
    f.create_dataset('markers/particle/id', data=ids, dtype='i8').attrs['unit'] = '1';

    f.close()

def write_hdf5_guidingcenters(fn, ids, mass, charge, r, phi, z, energy, pitch, theta, weight, time):
    """Write guiding center marker input in hdf5 file.

    Keyword arguments:
    fn     -- path to hdf5 file
    ids    -- unique identifier for each marker (positive integer)
    charge -- charge (e)
    mass   -- mass (amu)
    r      -- guiding center R coordinate (m)
    phi    -- guiding center phi coordinate (deg)
    z      -- guiding center z coordinate (m)
    energy -- guiding center energy (eV)
    pitch  -- guiding center pitch (v_para/v_tot)
    theta  -- guiding center gyroangle (rad)
    weight -- guiding center weight (markers/s)
    time   -- guiding center initial time (s)

    All guiding center inputs must be in numpy array format.
    """
    f = h5py.File(fn, "a")
    if not "/markers" in f:
        f.create_group('markers')

    if  "/markers/guiding_center" in f:
        del f["/markers/guiding_center"]


    f.create_group('markers/guiding_center')

    f['markers'].attrs['n_guiding_center'] = np.size(ids)

    if not 'n_particle' in f['markers'].attrs:
        f['markers'].attrs['n_particle'] = 0

    if not 'n_field_line' in f['markers'].attrs:
        f['markers'].attrs['n_field_line'] = 0

    f.create_dataset('markers/guiding_center/r', data=r, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset('markers/guiding_center/phi', data=phi, dtype='f8').attrs['unit'] = 'deg';
    f.create_dataset('markers/guiding_center/z', data=z, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset('markers/guiding_center/energy', data=energy, dtype='f8').attrs['unit'] = 'eV';
    f.create_dataset('markers/guiding_center/pitch', data=pitch, dtype='f8').attrs['unit'] = '1';
    f.create_dataset('markers/guiding_center/theta', data=theta, dtype='f8').attrs['unit'] = 'rad';
    f.create_dataset('markers/guiding_center/mass', data=mass, dtype='f8').attrs['unit'] = 'amu';
    f.create_dataset('markers/guiding_center/charge', data=charge, dtype='i4').attrs['unit'] = 'e';
    f.create_dataset('markers/guiding_center/weight', data=weight, dtype='f8').attrs['unit'] = 'markers/s';
    f.create_dataset('markers/guiding_center/time', data=time, dtype='f8').attrs['unit'] = 's';
    f.create_dataset('markers/guiding_center/id', data=ids, dtype='i8').attrs['unit'] = '1';

    f.close()

def write_hdf5_fieldlines(fn, ids, r, phi, z, pitch, weight, time):
    """Write magnetic field line marker input in hdf5 file.

    Keyword arguments:
    fn     -- path to hdf5 file
    ids    -- unique identifier for each marker (positive integer)
    r      -- magnetic field line R coordinate (m)
    phi    -- magnetic field line phi coordinate (deg)
    z      -- magnetic field line z coordinate (m)
    pitch  -- magnetic field line pitch whose sign defines 
              the direction field line is traced (real)
    weight -- magnetic field line weight (markers/s)
    time   -- magnetic field line initial time (s)

    All magnetic field line inputs must be in numpy array format.
    """
    f = h5py.File(fn, "a")
    if not "/markers" in f:
        f.create_group('markers')

    if  "/markers/field_line" in f:
        del f["/markers/field_line"]


    f.create_group('markers/field_line')

    f['markers'].attrs['n_field_line'] = np.size(ids)

    if not 'n_particle' in f['markers'].attrs:
        f['markers'].attrs['n_particle'] = 0

    if not 'n_guiding_center' in f['markers'].attrs:
        f['markers'].attrs['n_guiding_center'] = 0

    f.create_dataset('markers/field_line/r', data=r, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset('markers/field_line/phi', data=phi, dtype='f8').attrs['unit'] = 'deg';
    f.create_dataset('markers/field_line/z', data=z, dtype='f8').attrs['unit'] = 'm';
    f.create_dataset('markers/field_line/pitch', data=pitch, dtype='f8').attrs['unit'] = '1';
    f.create_dataset('markers/field_line/weight', data=weight, dtype='f8').attrs['unit'] = 'markers/s';
    f.create_dataset('markers/field_line/time', data=time, dtype='f8').attrs['unit'] = 's';
    f.create_dataset('markers/field_line/id', data=ids, dtype='i8').attrs['unit'] = '1';

    f.close()
