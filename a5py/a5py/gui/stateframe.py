"""
Contains definition of StateFrame class.

File: stateframe.py
"""
import copy
import tkinter
import tkinter.ttk as ttk

from mpl_toolkits.mplot3d import Axes3D

from .plotframe import PlotFrame
import a5py.marker.endcond as endcond

class StateFrame(PlotFrame):
    """
    A frame for plotting ini/endstate data.

    Side bar controls which quantities are plotted, are they on logarithmic axis
    and are they filtered by endconditions and are the plot axes scaled. Plot
    is replotted each time a control is changed. Axes are set automatically to
    3D if z coordinate is not None.
    """

    def __init__(self, gui, inistate, endstate=None):
        """
        Initialize and show default plot.
        """

        super().__init__(gui)
        self._inistate = inistate
        self._endstate = endstate
        self._states = ["inistate"]
        if endstate is not None:
            self._states.append("endstate")

        # Variables for plotting scatter plots.
        self._xchoice       = tkinter.StringVar(self)
        self._ychoice       = tkinter.StringVar(self)
        self._zchoice       = tkinter.StringVar(self)
        self._cchoice       = tkinter.StringVar(self)

        self._xlogchoice    = tkinter.IntVar(self)
        self._ylogchoice    = tkinter.IntVar(self)
        self._zlogchoice    = tkinter.IntVar(self)
        self._clogchoice    = tkinter.IntVar(self)

        self._equalchoice   = tkinter.IntVar(self)

        # Variables for plotting histograms.
        self._binxchoice    = tkinter.StringVar(self)
        self._binxlogchoice = tkinter.IntVar(self)
        self._binxminchoice = tkinter.IntVar(self)
        self._binxmaxchoice = tkinter.IntVar(self)
        self._nbinxchoice   = tkinter.IntVar(self)
        self._binychoice    = tkinter.StringVar(self)
        self._binylogchoice = tkinter.IntVar(self)
        self._binyminchoice = tkinter.IntVar(self)
        self._binymaxchoice = tkinter.IntVar(self)
        self._nbinychoice   = tkinter.IntVar(self)
        self._binzlogchoice = tkinter.IntVar(self)

        # These are needed by both.
        self._endcondchoice = tkinter.StringVar(self)
        self._statechoice   = tkinter.StringVar(self)
        self._plottype      = tkinter.StringVar(self)

        # List of all possible coordinates.
        self._coords = ["R", "phimod", "z", "time", "energy", "pitch", "vnorm",
                        "Bnorm", "vR", "vphi", "vz", "BR", "Bphi", "Bz", "mu",
                        "vpar", "charge", "id", "x", "y", "phi", "rho",
                        "polmod", "pol", "None"]

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
        self._binxchoice.set("R")
        self._binxlogchoice.set(0)
        self._binxminchoice.set(0)
        self._binxmaxchoice.set(20)
        self._nbinxchoice.set(10)
        self._binychoice.set("None")
        self._binylogchoice.set(0)
        self._binyminchoice.set(1)
        self._binymaxchoice.set(2)
        self._nbinychoice.set(10)
        self._binzlogchoice.set(0)
        self._endcondchoice.set("all")
        self._statechoice.set("inistate")
        self._plottype.set("scatter")

        # Update scatter plot each time a variable has changed.
        # (Histogram is updated when a button is pressed.)
        self._xchoice.trace(        'w', self._plot)
        self._ychoice.trace(        'w', self._plot)
        self._zchoice.trace(        'w', self._plot)
        self._cchoice.trace(        'w', self._plot)
        self._xlogchoice.trace(     'w', self._plot)
        self._ylogchoice.trace(     'w', self._plot)
        self._zlogchoice.trace(     'w', self._plot)
        self._clogchoice.trace(     'w', self._plot)
        self._equalchoice.trace(    'w', self._plot)
        self._endcondchoice.trace(  'w', self._plot)
        self._statechoice.trace(    'w', self._plot)

        # Make top panel to choose what kind of plot is shown.
        top = self.get_toppanel()
        buttona = tkinter.Radiobutton(top, text="Scatter",
                                      variable=self._plottype,
                                      value="scatter",
                                      command=self._showscatterpanel)
        buttonb = tkinter.Radiobutton(top, text="Histogram",
                                      variable=self._plottype,
                                      value="histogram",
                                      command=self._showhistogrampanel)
        buttona.pack(side="left")
        buttonb.pack(side="left")
        self._showscatterpanel()


    def _showscatterpanel(self):
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


    def _showhistogrampanel(self):
        pass


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        if self._statechoice.get() == "inistate":
            state = self._inistate
        else:
            state = self._endstate

        if self._plottype.get() == "scatter":
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
            if ccoord == "None":
                ccoord = None

            if zcoord == "None":
                axes = fig.add_subplot(1,1,1)
                state.scatter(xcoord, ycoord, c=ccoord, equal=equal,
                              log=log, endcond=endcond, axes=axes)
            else:
                # This would work for newer versions of matplotlib
                #ax = fig.add_subplot(1,1,1, projection="3d")
                axes = Axes3D(fig)
                state.scatter(xcoord, ycoord, zcoord, ccoord,
                              equal=equal, log=log,
                              endcond=endcond, axes=axes)

        else:
            pass

        self.draw()
