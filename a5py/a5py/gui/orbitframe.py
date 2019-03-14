"""
Contains definition of OrbitFrame class.

File: orbitframe.py
"""
import copy
import tkinter
import tkinter.ttk as ttk

import numpy as np

from mpl_toolkits.mplot3d import Axes3D

from .plotframe import PlotFrame
import a5py.marker.endcond as endcond

class OrbitFrame(PlotFrame):
    """
    A frame for plotting orbit data.

    Side bar controls which quantities are plotted, are they on logarithmic axis
    and are they filtered by endconditions and are the plot axes scaled. Plot
    is replotted each time a control is changed. Axes are set automatically to
    3D if z coordinate is not None.
    """

    def __init__(self, gui, orbits):
        """
        Initialize and show default plot.
        """
        super().__init__(gui)
        self._orbits = orbits

        # Variables for plotting orbits.
        self._xchoice       = tkinter.StringVar(self)
        self._ychoice       = tkinter.StringVar(self)
        self._zchoice       = tkinter.StringVar(self)
        self._cchoice       = tkinter.StringVar(self)

        self._xlogchoice    = tkinter.IntVar(self)
        self._ylogchoice    = tkinter.IntVar(self)
        self._zlogchoice    = tkinter.IntVar(self)
        self._clogchoice    = tkinter.IntVar(self)

        self._equalchoice   = tkinter.IntVar(self)

        # Variables for plotting Poincares.
        self._idtimechoice    = tkinter.IntVar(self)
        self._pncrcoordchoice = tkinter.StringVar(self)
        self._pncridchoice    = tkinter.IntVar(self)

        # These are needed by both.
        self._endcondchoice = tkinter.StringVar(self)
        self._plottype      = tkinter.StringVar(self)

        # List of all possible coordinates. Those not applicable to this data
        # will be removed.
        self._coords = ["R", "phimod", "z", "time", "energy", "pitch", "vnorm",
                        "Bnorm", "vR", "vphi", "vz", "BR", "Bphi", "Bz", "mu",
                        "vpar", "charge", "id", "x", "y", "phi", "rho",
                        "polmod", "pol", "None"]
        clist = copy.copy(self._coords)
        clist.remove("None")
        for c in clist:
            try:
                orbits[c]
            except (ValueError, AssertionError):
                self._coords.remove(c)

        # Check if this data contain Poincare data (and store Poincare ids).
        try:
            self._pncrids    = np.sort(np.unique(orbits["pncrid"])).tolist()
            self._pncrcoords = ["R-z", "rho-phi", "rho-theta", "R-phi"]
            haspoincare = True
        except (ValueError, AssertionError):
            haspoincare = False

        # All end conditions.
        self._endconds = ["all"] + list(endcond.endconds.keys())

        # Set default values for the variables.
        self._xchoice.set("R")
        self._ychoice.set("z")
        self._zchoice.set("None")
        self._cchoice.set("None")
        self._xlogchoice.set(0)
        self._ylogchoice.set(0)
        self._zlogchoice.set(0)
        self._clogchoice.set(0)
        self._equalchoice.set(1)
        self._idtimechoice.set(0)
        self._pncrcoordchoice.set("R-z")
        self._pncridchoice.set(0)
        self._endcondchoice.set("all")

        # Update plot each time a variable has changed.
        self._xchoice.trace(        'w', self._plot)
        self._ychoice.trace(        'w', self._plot)
        self._zchoice.trace(        'w', self._plot)
        self._cchoice.trace(        'w', self._plot)
        self._xlogchoice.trace(     'w', self._plot)
        self._ylogchoice.trace(     'w', self._plot)
        self._zlogchoice.trace(     'w', self._plot)
        self._clogchoice.trace(     'w', self._plot)
        self._equalchoice.trace(    'w', self._plot)
        self._idtimechoice.trace(   'w', self._plot)
        self._pncrcoordchoice.trace('w', self._plot)
        self._pncridchoice.trace(   'w', self._plot)
        self._endcondchoice.trace(  'w', self._plot)

        if haspoincare:
            self._plottype.set("poincare")
            top = self.get_toppanel()
            buttona = tkinter.Radiobutton(top, text="Orbits",
                                          variable=self._plottype,
                                          value="orbit",
                                          command=self._showorbitpanel)
            buttonb = tkinter.Radiobutton(top, text="Poincare",
                                          variable=self._plottype,
                                          value="poincare",
                                          command=self._showpoincarepanel)
            buttona.pack(side="left")
            buttonb.pack(side="left")
            self._showpoincarepanel()
        else:
            self._plottype.set("orbit")
            self._showorbitpanel()


    def _showorbitpanel(self):

        panel = self.get_sidepanel()

        xpanel = tkinter.Frame(panel)
        xlabel = tkinter.Label(xpanel, text="x: ")
        xinput = ttk.Combobox(xpanel, width=6, textvariable=self._xchoice)
        xltick = tkinter.Checkbutton(xpanel, text="log10", onvalue=1,
                                     offvalue=0, variable=self._xlogchoice,
                                     height=1, width=5)

        ypanel = tkinter.Frame(panel)
        ylabel = tkinter.Label(ypanel, text="y: ")
        yinput = ttk.Combobox(ypanel, width=6, textvariable=self._ychoice)
        yltick = tkinter.Checkbutton(ypanel, text="log10", onvalue=1,
                                     offvalue=0, variable=self._ylogchoice,
                                     height=1, width=5)

        zpanel = tkinter.Frame(panel)
        zlabel = tkinter.Label(zpanel, text="z: ")
        zinput = ttk.Combobox(zpanel, width=6, textvariable=self._zchoice)
        zltick = tkinter.Checkbutton(zpanel, text="log10", onvalue=1,
                                     offvalue=0, variable=self._zlogchoice,
                                     height=1, width=5)

        cpanel = tkinter.Frame(panel)
        clabel = tkinter.Label(cpanel, text="c: ")
        cinput = ttk.Combobox(cpanel, width=6, textvariable=self._cchoice)
        cltick = tkinter.Checkbutton(cpanel, text="log10", onvalue=1,
                                     offvalue=0, variable=self._clogchoice,
                                     height=1, width=5)

        tickequal = tkinter.Checkbutton(panel, text="Axis equal",
                                        variable=self._equalchoice,
                                        onvalue=1, offvalue=0,
                                        height=1, width=8)

        ecpanel = tkinter.Frame(panel)
        eclabel = tkinter.Label(ecpanel, text="Endcond: ")
        endcondinput = ttk.Combobox(ecpanel, width=6,
                                    textvariable=self._endcondchoice)

        xinput["values"] = self._coords
        yinput["values"] = self._coords
        zinput["values"] = self._coords
        cinput["values"] = self._coords
        endcondinput["values"] = self._endconds

        xlabel.pack(side="left")
        xinput.pack(side="left")
        xltick.pack(side="left")
        ylabel.pack(side="left")
        yinput.pack(side="left")
        yltick.pack(side="left")
        zlabel.pack(side="left")
        zinput.pack(side="left")
        zltick.pack(side="left")
        clabel.pack(side="left")
        cinput.pack(side="left")
        cltick.pack(side="left")

        eclabel.pack(side="left")
        endcondinput.pack(side="left")

        xpanel.pack()
        ypanel.pack()
        zpanel.pack()
        cpanel.pack()
        tickequal.pack()
        ecpanel.pack()

        self._plot()


    def _showpoincarepanel(self):
        panel = self.get_sidepanel()

        coordlabel = tkinter.Label(panel, text="Coordinates: ")
        coordinput = ttk.Combobox(panel, width=6,
                                  textvariable=self._pncrcoordchoice)
        idtimetick = tkinter.Checkbutton(panel, text="Id/Time", onvalue=1,
                                      offvalue=0, variable=self._idtimechoice,
                                      height=1, width=10)

        eclabel = tkinter.Label(panel, text="Endcond: ")
        endcondinput = ttk.Combobox(panel, width=6,
                                    textvariable=self._endcondchoice)

        pncrlabel = tkinter.Label(panel, text="Poincare: ")
        pncrinput = ttk.Combobox(panel, width=6,
                                 textvariable=self._pncridchoice)

        coordinput["values"]   = self._pncrcoords
        pncrinput["values"]    = self._pncrids
        endcondinput["values"] = self._endconds

        coordlabel.grid(row=0, column=0)
        coordinput.grid(row=0, column=1)
        idtimetick.grid(row=0, column=2)

        eclabel.grid(row=1, column=0)
        endcondinput.grid(row=1, column=1)

        pncrlabel.grid(row=2, column=0)
        pncrinput.grid(row=2, column=1)

        self._plot()


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        if self._plottype.get() == "orbit":
            xcoord  = self._xchoice.get()
            ycoord  = self._ychoice.get()
            zcoord  = self._zchoice.get()
            ccoord  = self._cchoice.get()
            equal   = self._equalchoice.get() == 1
            xlog    = self._xlogchoice.get()
            ylog    = self._ylogchoice.get()
            zlog    = self._zlogchoice.get()
            clog    = self._clogchoice.get()
            endcond = self._endcondchoice.get()

            if endcond == "all":
                endcond = None

            log = (xlog==1, ylog==1, zlog==1, clog==1)

            if xcoord == "None":
                xcoord = None
            if ycoord == "None":
                ycoord = None

            if zcoord == "None":
                axes = fig.add_subplot(1,1,1)

                if ccoord == "None":
                    self._orbits.plot(xcoord, ycoord, equal=equal, log=log,
                                      endcond=endcond, axes=axes)
                else:
                    self._orbits.scatter(xcoord, ycoord, c=ccoord, equal=equal,
                                         log=log, endcond=endcond, axes=axes)
            else:
                # This would work for newer versions of matplotlib
                #ax = fig.add_subplot(1,1,1, projection="3d")
                axes = Axes3D(fig)
                if ccoord == "None":
                    self._orbits.plot(xcoord, ycoord, zcoord, equal=equal,
                                      log=log, endcond=endcond, axes=axes)
                else:
                    self._orbits.scatter(xcoord, ycoord, zcoord, ccoord,
                                         equal=equal, log=log,
                                         endcond=endcond, axes=axes)

        else:
            # Coordinates are in form "x-y".
            xy = self._pncrcoordchoice.get()
            xcoord = xy.split("-")[0]
            ycoord = xy.split("-")[1]

            if xcoord in ["phi", "pol"]:
                xcoord += "mod"
            if ycoord in ["phi", "pol"]:
                ycoord += "mod"

            equal = False
            if xy == "R-z":
                equal = True

            endcond = self._endcondchoice.get()
            if endcond == "all":
                endcond = None

            pncrid = self._pncridchoice.get()

            axes = fig.add_subplot(1,1,1)

            if self._idtimechoice.get():
                self._orbits.poincare(xcoord, ycoord, pncrid, equal=equal,
                                      endcond=endcond, axes=axes)
            else:
                self._orbits.poincare(xcoord, ycoord, "time", pncrid,
                                      log=(0, 0, 0, 1), equal=equal,
                                      endcond=endcond, axes=axes)

        self.draw()
