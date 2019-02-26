"""
Helper routines and attributes for test modules.

File: testascot/helpers.py
"""
import subprocess
import numpy as np

from time import sleep
from time import perf_counter as timer

import a5py.ascot5io.ascot5      as ascot5
import a5py.ascot5io.ascot5tools as tools

## Name of the test HDF5 file
testfn = "test_ascot.h5"

## Name of the ascot5 binary
testbin = "ascot5_main"

def clean_opt(odict):
    """
    Turn all tunable features off in given options dictionary.

    Args:
        odict: dict Options dictionary
    """
    odict["RECORD_GO_AS_GC"] = np.array([0],dtype='i4')
    for o in odict:
        if o.startswith("ENABLE"):
            odict[o] = np.array([0],dtype='i4')
        if o.startswith("ENDCOND"):
            odict[o] = np.array([0],dtype='i4')
        if o.startswith("DISABLE"):
            odict[o] = np.array([0],dtype='i4')

def set_correct_input(parent, test):
    """
    Set that input field active that contains test name in its description.

    Args:
        parent:  str Input parent
        test:   str Name of the tests
    """
    a5 = ascot5.Ascot(testfn)
    typ = a5[parent][test].get_type()
    qid = a5[parent][test].get_qid()
    group = typ + "-" + qid
    tools.call_ascot5file(testfn, "set_active", group)

def set_and_run(test):
    """
    Set correct inputs active and carry out the test simulation.

    Args:
        testfn: str Name of the HDF5 file with tests
        test:   str Name of the test
    """
    set_correct_input("bfield",  test)
    set_correct_input("efield",  test)
    set_correct_input("marker",  test)
    set_correct_input("plasma",  test)
    set_correct_input("neutral", test)
    set_correct_input("wall",    test)
    set_correct_input("options", test)

    sleep(1.01) # Sleep for one second so that each run gets unique QID

    frm   = lambda x: "%.3f s" % x
    start = timer()
    subprocess.call(["./"+testbin, "--in="+testfn[:-3], "--d="+test],
                    stdout=subprocess.DEVNULL)
    print("Completed test " + test + " in " + frm(timer() - start))

def write_N0_3D_dummy(h5fn, desc):
    N0Rmin = 0
    N0Rmax = 100
    N0nR   = 2
    N0zmin = -100
    N0zmax = 100
    N0nz   = 2
    N0pmin = 0
    N0pmax = 2*np.pi
    N0np   = 2
    N0spec = 1
    N0anum = 1
    N0znum = 1
    N0dens = np.array([ [ [ [0,0] , [0,0] ], [ [0,0] , [0,0] ] ] ])
    N0_3D.write_hdf5(h5fn,
                     N0Rmin, N0Rmax, N0nR,
                     N0zmin, N0zmax, N0nz,
                     N0pmin, N0pmax, N0np,
                     N0spec, N0anum, N0znum, N0dens,
                     desc=desc)
