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

    if "/options" in f:
        del f["options"]
    if "/bfield" in f:
        del f["bfield"]
    if "/efield" in f:
        del f["efield"]
    if "/plasma" in f:
        del f["plasma"]
    f.close()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        clean(sys.argv[1])
    else:
        clean("ascot.h5")
