import sys
import h5py

def clean(fn):
    f = h5py.File(fn, "a")
    if  "inistate" in f:
        del f["inistate"]
    if  "endstate" in f:
        del f["endstate"]
    if  "orbits" in f:
        del f["orbits"]
    if  "dist" in f:
        del f["dist"]
    f.close()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        clean(sys.argv[1])
    else:
        clean("ascot.h5")
