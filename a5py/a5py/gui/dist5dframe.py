"""
Contains definition of StateFrame class.

File: dist5dframe.py
"""
import tkinter
import copy
import tkinter.ttk as ttk

from .plotframe import PlotFrame


class Dist5DFrame(PlotFrame):
    """
    An opening frame where other frames can be accessed and HDF5 file modified.
    """

    def __init__(self, gui, dist):
        """
        Initialize
        """
        super().__init__(gui)
        self._dist = dist

        self._xchoice   = tkinter.StringVar(self)
        self._ychoice   = tkinter.StringVar(self)
        self._equalaxis = tkinter.IntVar(self)

        self._xchoice.trace('w', self._set_xcoord)
        self._ychoice.trace('w', self._set_ycoord)

        top = self.get_toppanel()

        self._plottype = tkinter.StringVar()
        buttona = tkinter.Radiobutton(top, text="1D (vpa, vpe)",
                                      variable=self._plottype, value="1d",
                                      command=self._update_sidepanel)
        buttonb = tkinter.Radiobutton(top, text="2D (vpa, vpe)",
                                      variable=self._plottype, value="2d",
                                      command=self._update_sidepanel)
        buttonc = tkinter.Radiobutton(top, text="1D (E, xi)",
                                      variable=self._plottype, value="1dExi",
                                      command=self._update_sidepanel)
        buttond = tkinter.Radiobutton(top, text="2D (E, xi)",
                                      variable=self._plottype, value="2dExi",
                                      command=self._update_sidepanel)

        buttona.pack(side="left")
        buttonb.pack(side="left")
        buttonc.pack(side="left")
        buttond.pack(side="left")

        self._plottype.set("2d")
        self._update_sidepanel()


    def _update_sidepanel(self):
        val = self._plottype.get()

        if val == "1d":
            self._show_1dparams()
        elif val == "2d":
            self._show_2dparams()

        self._plot()


    def _show_1dparams(self):
        data = self._dist.get_dist()

        coords = data["abscissae"]
        self._xcoord = coords[0]
        self._ycoord = None

        self._coords = copy.copy(coords)

        panel = self.get_sidepanel()

        plotbutton = tkinter.Button(panel, text="Plot", command=self._plot)

        xinput = ttk.Combobox(panel, width=6, textvariable=self._xchoice)

        xinput.bind("<<ComboboxSelected>>", self._set_xcoord)
        xinput["values"] = self._coords

        xinput.pack()
        plotbutton.pack()

        self._plotbutton = plotbutton

        self._xchoice.set(coords[0])
        self._ychoice.set(coords[0])

        del data

        self._plot()


    def _show_2dparams(self):
        data = self._dist.get_dist()

        coords = data["abscissae"]
        self._xcoord = coords[0]
        self._ycoord = coords[1]

        self._coords = copy.copy(coords)

        panel = self.get_sidepanel()

        plotbutton = tkinter.Button(panel, text="Plot", command=self._plot)

        panel1 = tkinter.Frame(panel)
        xinput = ttk.Combobox(panel1, width=6, textvariable=self._xchoice)
        yinput = ttk.Combobox(panel1, width=6, textvariable=self._ychoice)

        tickequal = tkinter.Checkbutton(panel, text="Axis equal",
                                        variable=self._equalaxis,
                                        onvalue=1, offvalue=0,
                                        height=5, width=20)

        xinput.bind("<<ComboboxSelected>>", self._set_xcoord)
        xinput["values"] = self._coords

        yinput.bind("<<ComboboxSelected>>", self._set_ycoord)
        yinput["values"] = self._coords

        xinput.pack(side="left")
        yinput.pack(side="left")

        panel1.pack()
        tickequal.pack()
        plotbutton.pack()

        self._plotbutton = plotbutton

        self._xchoice.set(coords[0])
        self._ychoice.set(coords[2])
        self._equalaxis.set(1)

        del data

        self._plot()

    def _set_xcoord(self, *args):
        self._xcoord = self._xchoice.get()

        if (self._plottype.get() == "2d") and (self._xcoord == self._ycoord):
            self._plotbutton.config(state="disable")
        else:
            self._plotbutton.config(state="normal")

    def _set_ycoord(self, *args):
        self._ycoord = self._ychoice.get()

        if (self._plottype.get() == "2d") and (self._xcoord == self._ycoord):
            self._plotbutton.config(state="disable")
        else:
            self._plotbutton.config(state="normal")

    def _plot(self):
        fig = self.get_fig()
        ax = fig.add_subplot(1,1,1)

        if(self._plottype.get() == "1d"):
            self._dist.plot_dist(self._xcoord, axes=ax)

        elif(self._plottype.get() == "2d"):
            equal = self._equalaxis.get() == 1
            self._dist.plot_dist(self._xcoord, self._ycoord, equal=equal,
                                 axes=ax)

        self.draw()
