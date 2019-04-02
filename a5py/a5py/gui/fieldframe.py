"""
Contains definition of FieldFrame class.

File: fieldframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from .plotframe import PlotFrame
from .components import NumEntry

class FieldFrame(PlotFrame):
    """
    A frame for plotting any field-like data (scalara quantities on R,phi,z,t).
    """

    def __init__(self, gui, ascotpy, quantities):
        """
        Initialize and show default plot.
        """
        super().__init__(gui)
        self._gui    = gui
        self.ascotpy = ascotpy

        self._qchoice = tkinter.StringVar(self)

        # Set default values for the variables.
        self._qchoice.set(quantities[0])

        self.init2dpanel(quantities)


    def init2dpanel(self, quantities):
        panel = self.get_sidepanel()

        self.xmin_entry = NumEntry(panel, labeltext="R_min [m]:", defval=0.1)
        self.xmax_entry = NumEntry(panel, labeltext="R_max [m]:", defval=10)
        self.xnum_entry = NumEntry(panel, labeltext="R_num :",    defval=50,
                                   isint=True)

        self.ymin_entry = NumEntry(panel, labeltext="z_min [m]:", defval=-8)
        self.ymax_entry = NumEntry(panel, labeltext="z_max [m]:", defval=8)
        self.ynum_entry = NumEntry(panel, labeltext="z_num :",    defval=50,
                                   isint=True)

        self.xmin_entry.grid(row=0, column=0, sticky="W")
        self.xmax_entry.grid(row=1, column=0, sticky="W")
        self.xnum_entry.grid(row=2, column=0, sticky="W")

        self.ymin_entry.grid(row=0, column=1, sticky="W")
        self.ymax_entry.grid(row=1, column=1, sticky="W")
        self.ynum_entry.grid(row=2, column=1, sticky="W")

        self.phi_entry = NumEntry(panel, labeltext="phi [deg]:", defval=0)
        self.phi_entry.grid(row=3, column=0)

        self.time_entry = NumEntry(panel, labeltext="time [s]:", defval=0)
        self.time_entry.grid(row=4, column=0)

        qinput = ttk.Combobox(panel, width=10, textvariable=self._qchoice)
        qinput["values"] = quantities

        qinput.grid(row=6, column=1, sticky="WE")
        plotbutton = tkinter.Button(panel, text="Plot", command=self._plot)
        plotbutton.grid(row=8, column=1, sticky="WE")


    def _backtoindex(self):
        self.ascotpy.free(bfield=True, efield=True, plasma=True, neutral=True,
                          wall=True)
        super()._backtoindex()


    def _plot(self, *args):
        """
        Read control states and plot.

        Implement this in subclasses.
        """
        pass
