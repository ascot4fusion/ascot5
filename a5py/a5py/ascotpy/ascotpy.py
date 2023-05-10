"""
Python interface to ascot library.

This module defines Ascotpy which is a class that inherits all other
classes found in ascotpy package. Thus having an instance of Ascotpy is
enough to gain access to all available tools.

File: ascotpy.py
"""
import numpy as np

from a5py.ascotpy.libbfield   import LibBfield
from a5py.ascotpy.libefield   import LibEfield
from a5py.ascotpy.libplasma   import LibPlasma
from a5py.ascotpy.libneutral  import LibNeutral
from a5py.ascotpy.libboozer   import LibBoozer
from a5py.ascotpy.libmhd      import LibMhd
from a5py.ascotpy.libsimulate import LibSimulate

import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class Ascotpy(LibBfield, LibEfield, LibPlasma, LibNeutral, LibBoozer, LibMhd,
              LibSimulate):
    """
    One class to rule them all.
    """


    def evaluate(self, R, phi, z, t, quantity, grid=False,
                 squeeze=[None, None, None, None]):
        """
        Evaluate input quantities at given coordinates.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].
            grid : boolean, optional <br>
                treat input coordinates as abscissae and return the evaluated
                quantities on a grid instead.

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        R   = np.asarray(R).ravel()
        phi = np.asarray(phi).ravel()
        z   = np.asarray(z).ravel()
        t   = np.asarray(t).ravel()

        if grid:
            arrsize = (R.size, phi.size, z.size, t.size)
            R, phi, z, t = np.meshgrid(R, phi, z, t)
            R   = R.ravel()
            phi = phi.ravel()
            z   = z.ravel()
            t   = t.ravel()
        else:
            # Not a grid so check that dimensions are correct (and make
            # single-valued vectors correct size)
            arrsize = np.amax(np.array([R.size, phi.size, z.size, t.size]))
            errmsg = "Input arrays have inconsistent sizes ({}, {}, {}, {})"
            assert (R.size == 1 or R.size == arrsize) and  \
                (phi.size == 1 or phi.size == arrsize) and \
                (z.size == 1 or z.size == arrsize) and     \
                (t.size == 1 or t.size == arrsize),        \
                errmsg.format(R.size, phi.size, z.size, t.size)

            if R.size == 1:
                R = R[0]*np.ones((arrsize,))
            if phi.size == 1:
                phi = phi[0]*np.ones((arrsize,))
            if z.size == 1:
                z = z[0]*np.ones((arrsize,))
            if t.size == 1:
                t = t[0]*np.ones((arrsize,))

        out = None
        if quantity in LibBfield.quantities:
            out = LibBfield.evaluate(self, R, phi, z, t, quantity)
        if quantity in LibEfield.quantities:
            out = LibEfield.evaluate(self, R, phi, z, t, quantity)
        if quantity in LibNeutral.quantities:
            out = LibNeutral.evaluate(self, R, phi, z, t, quantity)
        if quantity in LibBoozer.quantities:
            out = LibBoozer.evaluate(self, R, phi, z, t, quantity)
        if quantity in LibMhd.quantities:
            out = LibMhd.evaluate(self, R, phi, z, t, quantity)
        if self.plasma_initialized and quantity in self.get_plasmaquantities():
            out = LibPlasma.evaluate(self, R, phi, z, t, quantity)

        if grid:
            out = np.reshape(out, arrsize)
            idx = [slice(1), slice(1), slice(1), slice(1)]
            for i in range(len(squeeze)):
                if squeeze[i] is not None:
                    # This does not work yet
                    out = np.apply_along_axis(squeeze[i], i, out)
                else:
                    idx[i] = slice(None)

            out = out[idx[0], idx[1], idx[2], idx[3]]

        return out


    def plotRz(self, R, phi, z, t, quantity, axes=None, clim=[None, None],
               **kwargs):
        """
        Evaluate and plot input quantities on Rz plane at given coordinates.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinate where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinate where data is evaluated [s].
            axes :  <br>
                plot on these axes instead of creating a new figure.
            clim : array_like <br>
                tuple with minimum and maximum color values.

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """

        out = self.evaluate(R, phi, z, t, quantity, grid=True)
        out = np.transpose(out[:,0,:,0])

        newfig = axes is None
        if newfig:
            plt.figure()
            axes = plt.gca()

        out = np.ma.masked_invalid(out)
        if clim[0] is None:
            clim[0] = np.nanmin(out)
        if clim[1] is None:
            clim[1] = np.nanmax(out)

        mesh = axes.pcolormesh( R, z, out, vmin=clim[0], vmax=clim[1] )
        axes.patch.set(hatch='x', edgecolor=[0.9, 0.9, 0.9])

        plt.colorbar(mesh)
        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(R[0], R[-1])
        axes.set_ylim(z[0], z[-1])

        if newfig:
            plt.show(block=False)

        return axes
