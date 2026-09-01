import numpy as np

from a5py.engine.interpolate import evaluate
from a5py.plotting import openfigureifnoaxes, mesh2d

@openfigureifnoaxes()
def plot_poloidal(ascot, qnt, phi=0, t=0, nr=100, nz=200, axes=None):
    scale = ascot.data.bfield.active.rmajor
    r = np.linspace(scale/2, scale*1.5, nr)
    z = np.linspace(-scale/2, scale/2, nz)
    val = evaluate(r, phi, z, t, qnt, grid=True, bfield=ascot.data.bfield.active)
    val = np.squeeze(val)
    mesh2d(r, z, val, axes=axes)


def plot_topdown(ascot):
    pass