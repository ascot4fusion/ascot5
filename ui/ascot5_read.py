import numpy as np
import h5py
import ascot5_orbits

def ascot5_read(fn):
    f = h5py.File(fn, "r")

    mode = 2

    # Read the requested input if present
    out = {}
    if "/options" in f:
        out["options"] = {} #ascot5_options_read(f["options"], mode)
    else:
        out["options"] = {}

    if "/bfield" in f:
        out["bfield"] = {} #ascot5_bfield_read(f["bfield"], mode)
    else:
        out["bfield"] = {}

    if "/efield" in f:
        out["efield"] = {} #ascot5_efield_read(f["efield"], mode)
    else:
        out["efield"] = {}

    if "/wall" in f:
        out["wall"] = {} #ascot5_wall_read(f["wall"], mode)
    else:
        out["wall"] = {}

    if "/plasma" in f:
        out["plasma"] = {} #ascot5_plasma_read(f["plasma"], mode)
    else:
        out["plasma"] = {}

    if "/comments" in f:
        out["comments"] = {} #ascot5_comments_read(f["comments"], mode)
    else:
        out["comments"] = {}

    if "/markers" in f:
        out["markers"] = {} #ascot5_markers_read(f["markers"], mode)
    else:
        out["markers"] = {}

    if "/dists" in f:
        out["dists"] = {} #ascot5_dists_read(f["dists"], mode)
    else:
        out["dists"] = {}

    if "/orbits" in f:
        out["orbits"] = ascot5_orbits.read(f["orbits"], mode)
    else:
        out["orbits"] = {}

    if "/states" in f:
        out["states"] = {} #ascot5_states_read(f["states"], mode)
    else:
        out["states"] = {}

    f.close()
    return out
