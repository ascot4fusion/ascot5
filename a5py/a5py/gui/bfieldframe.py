"""
Contains definition of BfieldFrame class.

File: bfieldframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from a5py.ascotpy.libbfield import LibBfield

from .fieldframe import FieldFrame

class BfieldFrame(FieldFrame):
    """
    A frame for plotting magnetic field data.
    """

    def __init__(self, gui, ascotpy, walldata=None):
        """
        Initialize and show default plot.
        """
        ascotpy.init(bfield=True)
        super().__init__(gui, ascotpy, LibBfield.quantities)

        self._walldata = walldata

        self._plot()


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

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
