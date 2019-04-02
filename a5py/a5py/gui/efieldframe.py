"""
Contains definition of EfieldFrame class.

File: efieldframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from a5py.ascotpy.libefield import LibEfield

from .fieldframe import FieldFrame

class EfieldFrame(FieldFrame):
    """
    A frame for plotting electric field data.
    """

    def __init__(self, gui, ascotpy, walldata=None):
        """
        Initialize and show default plot.
        """
        ascotpy.init(bfield=True, efield=True)
        super().__init__(gui, ascotpy, LibEfield.quantities)

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

        axes = fig.add_subplot(1,1,1)
        self.ascotpy.plotRz(r, phi, z, time, self._qchoice.get(), axes)
        self.ascotpy.plotseparatrix(r, phi, z, time, axes)

        if self._walldata is not None:
            self._walldata.plotRz(axes=axes, phi=phi)

        self.draw()
