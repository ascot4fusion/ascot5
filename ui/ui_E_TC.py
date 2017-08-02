import h5py
import numpy as np
	
def write_hdf5(fn, Exyz):
    """Write trivial cartesian electric field input in hdf5 file.

    Keyword arguments:
    fn      -- path to hdf5 file
    Exyz    -- Electric field value in Cartesian coordinates (V/m)
    """

    # check that input is a valid array with three elements
    if Exyz.shape != (3,1) and Exyz.shape != (1,3) and Exyz.shape != (3,):
        raise Exception('Exyz has invalid format.')

    f = h5py.File(fn, "a")

    # clear existing groups
    if not "/efield" in f:
        o = f.create_group('efield')
    else:
        o = f["efield"]
        del o.attrs["type"]

    if  "/efield/E_TC" in f:
        del f["/efield/E_TC"]

    # write dataset
    o.attrs["type"] = np.string_("E_TC")
    o = o.create_group("E_TC")
    o.create_dataset("Exyz", (3,1), dtype='f8', data = Exyz)
    f.close()
