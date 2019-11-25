#!/usr/bin/env python3

"""
Test ASCOT5 physics.

This module contains the routines to carry out all physics tests. These routines
call the init, run, and check routines of each test module.

All tests are done using a single HDF5 file test_ascot.h5. It is user's
responsibility to remove that file before repeating the tests. Running
tests requires that ascot5_main binary is in the same folder as from which
this script is called.

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
from time import perf_counter as timer

import numpy as np

import a5py.ascot5io.ascot5      as ascot5
import a5py.ascot5io.ascot5tools as tools

import a5py.testascot.helpers as helpers
import test_elementary
import test_orbitfollowing
import test_gctransform
import test_coulombcollisions
import test_classicaltransport
import test_neoclassicaltransport

# You can specify here which tests to run

## Test test_elementary.py
dotest_elementary            = True
## Test test_orbitfollowing.py
dotest_orbitfollowing        = True
## Test test_gctransform.py
dotest_gctransform           = True
## Test test_coulombcollisions.py
dotest_coulombcollisions     = True
## Test test_classicaltransport.py
dotest_classicaltransport    = True
## Test test_neoclassicaltransport.py
dotest_neoclassicaltransport = True

def init():
    """
    Initialize all tests to test_ascot.h5
    """

    # Make sure no existing data interferes with new tests.
    subprocess.call(["rm", "-f", helpers.testfn])

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
    if not os.path.isfile(helpers.testfn):
        print(helpers.testfn + " is missing.")
        print("Aborting runs.")
        sys.exit()

    if not os.path.isfile(helpers.testbin):
        print(helpers.testbin + " is missing.")
        print("Aborting runs.")
        sys.exit()

    frm   = lambda x: "%.3f s" % x
    start = timer()
    inter = start

    if dotest_elementary:
        print("Running test_elementary.")
        test_elementary.run()
        print("Elapsed time is " + frm(timer() - inter))
        print("")
        inter = timer()

    if dotest_orbitfollowing:
        print("Running test_orbitfollowing.")
        test_orbitfollowing.run()
        print("Elapsed time is " + frm(timer() - inter))
        print("")
        inter = timer()

    if dotest_gctransform:
        print("Running test_gctransform.")
        test_gctransform.run()
        print("Elapsed time is " + frm(timer() - inter))
        print("")
        inter = timer()

    if dotest_coulombcollisions:
        print("Running test_coulombcollisions.")
        test_coulombcollisions.run()
        print("Elapsed time is " + frm(timer() - inter))
        print("")
        inter = timer()

    if dotest_classicaltransport:
        print("Running test_classicaltransport.")
        test_classicaltransport.run()
        print("Elapsed time is " + frm(timer() - inter))
        print("")
        inter = timer()

    if dotest_neoclassicaltransport:
        print("Running test_neoclassicaltransport.")
        test_neoclassicaltransport.run()
        print("Elapsed time is " + frm(timer() - inter))
        print("")
        inter = timer()

    print("Total elapsed time is " + frm(timer() - start))


def check():
    """
    Plot the results of these tests.
    """
    if not os.path.isfile(helpers.testfn):
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
