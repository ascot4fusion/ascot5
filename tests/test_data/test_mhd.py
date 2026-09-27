import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.engine.interpolate import evaluate


def test_mhdstationary_init():
    ascot = Ascot()

    rhogrid = np.linspace(0, 1, 10)
    toroidalnumber = [2, 5]
    poloidalnumber = [3, 6]
    phase = np.random.random((len(toroidalnumber),)) * unyt.rad
    frequency = np.random.random((len(toroidalnumber),)) * unyt.rad/unyt.s
    amplitude = np.random.random((len(toroidalnumber),))
    magneticprofile = np.random.random((rhogrid.size, len(toroidalnumber))) * unyt.m
    electricprofile = np.random.random((rhogrid.size, len(toroidalnumber))) * unyt.V
    obj = ascot.data.create_mhdstationary(
        rhogrid=rhogrid,
        toroidalnumber=toroidalnumber,
        poloidalnumber=poloidalnumber,
        phase=phase,
        frequency=frequency,
        amplitude=amplitude,
        magneticprofile=magneticprofile,
        electricprofile=electricprofile,
    )

    assert np.allclose(obj.rhogrid, rhogrid)
    assert np.allclose(obj.toroidalnumber, toroidalnumber)
    assert np.allclose(obj.poloidalnumber, poloidalnumber)
    assert np.allclose(obj.magneticprofile, magneticprofile)
    assert np.allclose(obj.electricprofile, electricprofile)
    assert np.allclose(obj.phase, phase)
    assert np.allclose(obj.frequency, frequency)
    assert np.allclose(obj.amplitude, amplitude)

def test_mhddynamic_init():
    ascot = Ascot()

    rhogrid = np.linspace(0, 1, 3)
    timegrid = np.linspace(0, 10, 4) * unyt.s
    toroidalnumber = [2, 5]
    poloidalnumber = [3, 6]
    magneticprofile = np.random.random((rhogrid.size, timegrid.size, len(toroidalnumber))) * unyt.m
    electricprofile = np.random.random((rhogrid.size, timegrid.size, len(toroidalnumber))) * unyt.V
    phase = np.random.random((len(toroidalnumber),)) * unyt.rad
    frequency = np.random.random((len(toroidalnumber),)) * unyt.rad/unyt.s
    amplitude = np.random.random((len(toroidalnumber),))
    obj = ascot.data.create_mhddynamic(
        rhogrid=rhogrid,
        timegrid=timegrid,
        toroidalnumber=toroidalnumber,
        poloidalnumber=poloidalnumber,
        magneticprofile=magneticprofile,
        electricprofile=electricprofile,
        phase=phase,
        frequency=frequency,
        amplitude=amplitude
    )

    assert np.allclose(obj.rhogrid, rhogrid)
    assert np.allclose(obj.toroidalnumber, toroidalnumber)
    assert np.allclose(obj.poloidalnumber, poloidalnumber)
    assert np.allclose(obj.magneticprofile, magneticprofile)
    assert np.allclose(obj.electricprofile, electricprofile)
    assert np.allclose(obj.phase, phase)
    assert np.allclose(obj.frequency, frequency)
    assert np.allclose(obj.amplitude, amplitude)
