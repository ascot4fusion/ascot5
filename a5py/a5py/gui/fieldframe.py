"""
Contains definition of FieldFrame class.

File: fieldframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from .plotframe import PlotFrame

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

        self._binxlogchoice = tkinter.IntVar(self)
        self._binxminchoice = tkinter.DoubleVar(self)
        self._binxmaxchoice = tkinter.DoubleVar(self)
        self._nbinxchoice   = tkinter.IntVar(self)
        self._binylogchoice = tkinter.IntVar(self)
        self._binyminchoice = tkinter.DoubleVar(self)
        self._binymaxchoice = tkinter.DoubleVar(self)
        self._nbinychoice   = tkinter.IntVar(self)
        self._phichoice     = tkinter.DoubleVar(self)

        self._qchoice = tkinter.StringVar(self)

        # Set default values for the variables.
        self._binxminchoice.set(1)
        self._binxmaxchoice.set(20)
        self._nbinxchoice.set(100)
        self._binylogchoice.set(0)
        self._binyminchoice.set(-10)
        self._binymaxchoice.set(10)
        self._nbinychoice.set(100)
        self._phichoice.set(0)
        self._qchoice.set(quantities[0])

        self.init2dpanel(quantities)


    def init2dpanel(self, quantities):
        panel = self.get_sidepanel()

        vcmd = (self.register(self._validate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        binpanel   = tkinter.Frame(panel)
        xminlabel  = tkinter.Label(binpanel, text="R min ")
        xmaxlabel  = tkinter.Label(binpanel, text="R max ")
        nxlabel    = tkinter.Label(binpanel, text="R nbin")
        minxentry  = tkinter.Entry(binpanel, validate = 'key',
                                   validatecommand = vcmd, width=6,
                                   textvariable=self._binxminchoice)
        maxxentry  = tkinter.Entry(binpanel, validate = 'key',
                                   validatecommand = vcmd, width=6,
                                   textvariable=self._binxmaxchoice)
        nxentry    = tkinter.Entry(binpanel, validate = 'key',
                                   validatecommand = vcmd, width=6,
                                   textvariable=self._nbinxchoice)

        yminlabel   = tkinter.Label(binpanel, text="z min ")
        ymaxlabel   = tkinter.Label(binpanel, text="z max ")
        nylabel     = tkinter.Label(binpanel, text="z nbin")
        minyentry   = tkinter.Entry(binpanel, validate = 'key',
                                    validatecommand = vcmd, width=6,
                                    textvariable=self._binyminchoice)
        maxyentry   = tkinter.Entry(binpanel, validate = 'key',
                                    validatecommand = vcmd, width=6,
                                    textvariable=self._binymaxchoice)
        nyentry     = tkinter.Entry(binpanel, validate = 'key',
                                    validatecommand = vcmd, width=6,
                                    textvariable=self._nbinychoice)
        philabel    = tkinter.Label(binpanel, text="phi ")
        phientry    = tkinter.Entry(binpanel, validate = 'key',
                                    validatecommand = vcmd, width=6,
                                    textvariable=self._phichoice)

        xminlabel.grid(row=0, column=0)
        minxentry.grid(row=0, column=1)
        xmaxlabel.grid(row=1, column=0)
        maxxentry.grid(row=1, column=1)
        nxlabel.grid(  row=2, column=0)
        nxentry.grid(  row=2, column=1)

        yminlabel.grid(row=0, column=2)
        minyentry.grid(row=0, column=3)
        ymaxlabel.grid(row=1, column=2)
        maxyentry.grid(row=1, column=3)
        nylabel.grid(  row=2, column=2)
        nyentry.grid(  row=2, column=3)
        philabel.grid( row=3, column=0)
        phientry.grid( row=3, column=1)

        binpanel.pack()

        qinput = ttk.Combobox(panel, width=10, textvariable=self._qchoice)
        qinput["values"] = quantities

        qinput.pack()
        tkinter.Button(panel, text="Plot", command=self._plot).pack()


    def _backtoindex(self):
        self.ascotpy.free(bfield=True)
        super()._backtoindex()


    def _validate(self, action, index, value_if_allowed,
                  prior_value, text, validation_type, trigger_type,
                  widget_name):
        if action == "1":
            if text in "0123456789.-+":
                try:
                    float(value_if_allowed)
                    return True
                except ValueError:
                    return False
            else:
                return False


    def _plot(self, *args):
        """
        Read control states and plot.

        Implement this in subclasses.
        """
        pass
