"""
Test ASCOT5 physics.

This module contains the routines to carry out all physics tests. These routines
call the init, run, and check routines of each test module. This module also
contains some helper routines required by the test modules.

All tests are done using a single HDF5 file test_ascot.h5. It is user's
responsibility to remove that file before repeating the tests. Running
tests requires that ascot5_main binary is in the same folder as from which
this script is called. The opt.py file is also required and it can be either on
a same folder (i.e. in test_ascot folder) or the parent folder.

To init, run and check all tests, call this script as
> python test_ascot.py

To only init (i.e. write input data to the HDF5 file), call
> python test_ascot.py init

To only run (i.e. carry out the ascot simulations), call
> python test_ascot.py run

To only check (i.e. plot the results), call
> python test_ascot.py check

File: test_ascot.py
"""

import sys
import subprocess
import os.path
from time import sleep

import numpy as np

import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.ascot5tools as tools

sys.path.insert(0, '../')
sys.path.insert(0, '.')
import test_elementary
import test_orbitfollowing
import test_gctransform
import test_coulombcollisions
import test_classicaltransport
import test_neoclassicaltransport

testfn = "test_ascot.h5"
testbin = "ascot5_main"

# You can specify here which tests to run

## Test test_elementary.py
dotest_elementary            = False
## Test test_orbitfollowing.py
dotest_orbitfollowing        = False
## Test test_gctransform.py
dotest_gctransform           = False
## Test test_coulombcollisions.py
dotest_coulombcollisions     = True
## Test test_classicaltransport.py
dotest_classicaltransport    = False
## Test test_neoclassicaltransport.py
dotest_neoclassicaltransport = False

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
    qid = a5[parent][test].get_qid()
    tools.set_active(testfn, qid)

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

    sleep(1.01)
    subprocess.call(["./"+testbin, "--in="+testfn[:-3], "--d="+test],
                    stdout=subprocess.DEVNULL)
    print("Completed test " + test)

def init():
    """
    Initialize all tests to test_ascot.h5
    """

    # Make sure no existing data interferes with new tests.
    subprocess.call(["rm", "-f", testfn])

    if dotest_elementary:
        print("Initializing test_elementary.")
        test_elementary.init()

    if dotest_orbitfollowing:
        print("Initializing test_orbitfollowing.")
        test_orbitfollowing.init()

    if dotest_gctransform:
        print("Initializing test_gctransform.")
        test_gctransform.init()

    if dotest_coulombcollisions:
        print("Initializing test_coulombcollisions.")
        test_coulombcollisions.init()

    if dotest_classicaltransport:
        print("Initializing test_classicaltransport.")
        test_classicaltransport.init()

    if dotest_neoclassicaltransport:
        print("Initializing test_neoclassicaltransport.")
        test_neoclassicaltransport.init()


def run():
    """
    Run each test

    This function goes through each test by first finding which inputs are
    required for that particular test (from their description) and using the
    qids to carry out the simulation. A new simulation is started when the
    previous one finishes.
    """
    if not os.path.isfile(testfn):
        print(testfn + " is missing.")
        print("Aborting runs.")
        sys.exit()

    if not os.path.isfile(testbin):
        print(testbin + " is missing.")
        print("Aborting runs.")
        sys.exit()

    if dotest_elementary:
        print("Running test_elementary.")
        test_elementary.run()

    if dotest_orbitfollowing:
        print("Running test_orbitfollowing.")
        test_orbitfollowing.run()

    if dotest_gctransform:
        print("Running test_gctransform.")
        test_gctransform.run()

    if dotest_coulombcollisions:
        print("Running test_coulombcollisions.")
        test_coulombcollisions.run()

    if dotest_classicaltransport:
        print("Running test_classicaltransport.")
        test_classicaltransport.run()

    if dotest_neoclassicaltransport:
        print("Running test_neoclassicaltransport.")
        test_neoclassicaltransport.run()


def check():
    """
    Plot the results of these tests.
    """
    if not os.path.isfile(testfn):
        print(testfn + " is missing.")
        print("Aborting checks.")
        sys.exit()

    if dotest_elementary:
        print("Checking results of test_elementary.")
        test_elementary.check()

    if dotest_orbitfollowing:
        print("Checking results of test_orbitfollowing.")
        test_orbitfollowing.check()

    if dotest_gctransform:
        print("Checking results of test_gctransform.")
        test_gctransform.check()

    if dotest_coulombcollisions:
        print("Checking results of test_coulombcollisions.")
        test_coulombcollisions.check()

    if dotest_classicaltransport:
        print("Checking results of test_classicaltransport.")
        test_classicaltransport.check()

    if dotest_neoclassicaltransport:
        print("Checking results of test_neoclassicaltransport.")
        test_neoclassicaltransport.check()

if __name__ == '__main__':

    if( len(sys.argv) == 1 ):
        print("Initializing tests.")
        init()
        print("Initialization complete.")
        print("")
        print("Running tests.")
        run()
        print("Runs complete.")
        print("")
        print("Checking test results.")
        check()
        print("Testing complete.")
        sys.exit()

    if(len(sys.argv) > 2):
        print("Too many arguments.")
        print("Only \"init\", \"run\" or \"check\" is accepted.")
        print("Aborting.")
        sys.exit()

    if(   sys.argv[1] == "init" ):
        print("Initializing tests.")
        init()
        print("Initialization complete.")
        sys.exit()

    elif( sys.argv[1] == "run" ):
        print("Running tests.")
        run()
        print("Runs complete.")
        sys.exit()

    elif( sys.argv[1] == "check" ):
        print("Checking test results.")
        check()
        print("Testing complete.")
        sys.exit()

    else:
        print("Too many arguments.")
        print("Only \"init\", \"run\" or \"check\" is accepted.")
        print("Aborting.")
        sys.exit()
