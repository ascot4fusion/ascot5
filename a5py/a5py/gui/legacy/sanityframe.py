"""
Contains definition of SanityFrame class.

File: sanityframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from .plotframe import PlotFrame
from .components import NumEntry

class SanityFrame(PlotFrame):
    """
    A frame for plotting fast sanity checks.
    """

    def __init__(self, gui, ascotpy, qid, marker):
        """
        Initialize and show plot.
        """
        super().__init__(gui)
        self._gui    = gui
        self.ascotpy = ascotpy
        self.ascotpy.init(bfield=qid)
        self.marker = marker

        self._plot()



    def _backtoindex(self):
        self.ascotpy.free(bfield=True)
        super()._backtoindex()


    def _plot(self, *args):
        axis = self.ascotpy.evaluate(1, 0, 0, 0, "axis")
        z0 = axis["axisz"]
        r0 = axis["axisr"]

        rmin = r0-r0/2
        rmax = r0+r0/2
        dphi = 10 * np.pi/180 # So that b and j quivers dont overlap

        r   = np.linspace(rmin, rmax, 10)
        phi = np.linspace(0, 360, 18, endpoint=False) * np.pi/180
        phi, r = np.meshgrid(phi, r)
        r   = r.ravel()
        phi = phi.ravel()

        br   = np.squeeze(self.ascotpy.evaluate(r, phi, z0, 0, "br"))
        bphi = np.squeeze(self.ascotpy.evaluate(r, phi, z0, 0, "bphi"))
        jr   = np.squeeze(self.ascotpy.evaluate(r, phi + dphi, z0, 0, "jr"))
        jphi = np.squeeze(self.ascotpy.evaluate(r, phi + dphi, z0, 0, "jphi"))

        x  = np.cos(phi) * r
        y  = np.sin(phi) * r
        bx = np.cos(phi) * br - np.sin(phi) * bphi
        by = np.sin(phi) * br + np.cos(phi) * bphi
        bnorm = np.sqrt(bx**2 + by**2)

        xj = np.cos(phi+dphi) * r
        yj = np.sin(phi+dphi) * r
        jx = np.cos(phi+dphi) * jr - np.sin(phi+dphi) * jphi
        jy = np.sin(phi+dphi) * jr + np.cos(phi+dphi) * jphi
        jnorm = np.sqrt(jx**2 + jy**2)

        fig = self.get_fig()

        axes = fig.add_subplot(1,1,1)
        axes.quiver(x,y, bx/bnorm, by/bnorm,
                    color="blue", scale=40)
        axes.quiver(xj,yj, jx/jnorm, jy/jnorm,
                    color="red", scale=40)

        axes.set_aspect("equal", adjustable="box")

        marker = self.marker.read()
        x  = np.cos(marker["phi"] * np.pi/180) * marker["r"]
        y  = np.sin(marker["phi"] * np.pi/180) * marker["r"]

        axes.scatter(x,y, s=1, c="black", zorder=-2)

        axes.legend(("Magnetic field", "Plasma current", "Markers"))

        self.draw()
