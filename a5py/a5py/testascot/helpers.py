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
import a5py.ascot5io.N0_3D       as N0_3D

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
    odict["RECORD_MODE"] = np.array([0],dtype='i4')
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
    set_correct_input("boozer",  test)
    set_correct_input("mhd",     test)
    
    sleep(1.01) # Sleep for one second so that each run gets unique QID

    frm   = lambda x: "%.3f s" % x
    start = timer()
    subprocess.call(["./"+testbin, "--in="+testfn[:-3], "--d="+test],
                    stdout=subprocess.DEVNULL)
    print("Completed test " + test + " in " + frm(timer() - start)) 
