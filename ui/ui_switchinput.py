import numpy as np
import h5py

if __name__ == "__main__":
    if len(sys.argv) > 1:
        switchinput(sys.argv[1], sys.argv[2], sys.argv[2])
    else:
        fn = 'acsot.h5'
        fieldname = 'bfield'
        fieldattr = 'B2D'
        switchinput(fn,fieldname, fieldattr)

def switchinput(fn,fieldname, fieldattr):
    """Switch between different types of input

    Example: ascot.h5 includes booth /bfield/B2D and
    /bfield/B3D. User can choose to use 2D field by calling

    switchinput('ascot.h5','bfield','B2D')

    and back to 3D as

    switchinput('ascot.h5','bfield','B3D')

    Keyword arguments:
    fn        -- path to hdf5 file
    fieldname -- input field whose type is affected
    fieldattr -- the new type
    """
    f = h5py.File(fn, "a")
    o = f[fieldname]
    del o.attrs["type"]
    o.attrs["type"] = np.string_(fieldattr)
    f.close()


    
