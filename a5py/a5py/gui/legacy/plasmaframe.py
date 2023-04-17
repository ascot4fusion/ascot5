"""
Contains definition of PlasmaFrame class.

File: plasmaframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from a5py.ascotpy.libplasma import LibPlasma

from .fieldframe import FieldFrame

class PlasmaFrame(FieldFrame):
    """
    A frame for plotting plasma data.
    """

    def __init__(self, gui, ascotpy, qid, bqid, walldata=None):
        """
        Initialize and show default plot.
        """
        ascotpy.init(bfield=bqid, plasma=qid)
        super().__init__(gui, ascotpy, ascotpy.get_plasmaquantities())

        self._walldata = walldata

        top = self.get_toppanel()
        self._plottype = tkinter.StringVar()
        buttona = tkinter.Radiobutton(top, text="1D", variable=self._plottype,
                                      value="1d", command=self._plot())
        buttonb = tkinter.Radiobutton(top, text="2D", variable=self._plottype,
                                      value="2d", command=self._plot())
        buttona.pack(side="left")
        buttonb.pack(side="left")
        self._plottype.set("2d")

        self._plot()


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        if self._plottype.get() == "1d":
            axis = self.ascotpy.evaluate(0, 0, 0, 0, "axis")
            r = np.linspace(axis["axisr"][0],10,1000)
            rhov = self.ascotpy.evaluate(r, 0, axis["axisz"][0], 0, "rho")
            quantity = self.ascotpy.evaluate(r[rhov <= 1], axis["axisz"][0], 0,
                                             0, self._qchoice.get())
            axes = fig.add_subplot(1,1,1)
            axes.plot(rhov[rhov <= 1], quantity)

        else:
            r = np.linspace( self.xmin_entry.getval() ,
                             self.xmax_entry.getval() ,
                             self.xnum_entry.getval() )

            z = np.linspace( self.ymin_entry.getval() ,
                             self.ymax_entry.getval() ,
                             self.ynum_entry.getval() )

            phi  = self.phi_entry.getval() * np.pi / 180
            time = self.time_entry.getval()

            clim = [None, None]
            if not self.cmin_entry.isempty():
                clim[0] = self.cmin_entry.getval()
            if not self.cmax_entry.isempty():
                clim[1] = self.cmax_entry.getval()

            axes = fig.add_subplot(1,1,1)
            self.ascotpy.plotRz(r, phi, z, time, self._qchoice.get(), axes=axes,
                                clim=clim)
            self.ascotpy.plotseparatrix(r, phi, z, time, axes)

            if self._walldata is not None:
                self._walldata.plotRz(axes=axes, phi=phi)

        self.draw()
