"""
Contains definition of OrbitFrame class.

File: orbitframe.py
"""
import copy
import tkinter
import tkinter.ttk as ttk

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

        self._xchoice       = tkinter.StringVar(self)
        self._ychoice       = tkinter.StringVar(self)
        self._zchoice       = tkinter.StringVar(self)

        self._equalchoice   = tkinter.IntVar(self)

        self._endcondchoice = tkinter.StringVar(self)

        self._xlogchoice    = tkinter.IntVar(self)
        self._ylogchoice    = tkinter.IntVar(self)
        self._zlogchoice    = tkinter.IntVar(self)

        # List of all possible coordinates. Those not applicable to this data
        # will be removed.
        coords = ["R", "phimod", "z", "time", "energy", "pitch", "vnorm",
                  "Bnorm", "vR", "vphi", "vz", "BR", "Bphi", "Bz", "mu", "vpar",
                  "charge", "id", "x", "y", "phi", "None"]
        clist = copy.copy(coords)
        clist.remove("None")
        for c in clist:
            try:
                orbits[c]
            except ValueError:
                coords.remove(c)

        endconds = ["all"] + list(endcond.endconds.keys())

        self._xchoice.set(coords[0])
        self._ychoice.set(coords[2])
        self._zchoice.set(coords[-1])
        self._xlogchoice.set(0)
        self._ylogchoice.set(0)
        self._zlogchoice.set(0)
        self._equalchoice.set(1)
        self._endcondchoice.set("all")

        self._xchoice.trace('w', self._plot)
        self._ychoice.trace('w', self._plot)
        self._zchoice.trace('w', self._plot)
        self._equalchoice.trace('w', self._plot)
        self._endcondchoice.trace('w', self._plot)
        self._xlogchoice.trace('w', self._plot)
        self._ylogchoice.trace('w', self._plot)
        self._zlogchoice.trace('w', self._plot)

        self._xcoord  = self._xchoice.get()
        self._ycoord  = self._ychoice.get()
        self._zcoord  = self._zchoice.get()
        self._equal   = self._equalchoice.get()
        self._xlog    = self._xlogchoice.get()
        self._ylog    = self._ylogchoice.get()
        self._zlog    = self._zlogchoice.get()
        self._endcond = self._endcondchoice.get()

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

        tickequal = tkinter.Checkbutton(panel, text="Axis equal",
                                        variable=self._equalchoice,
                                        onvalue=1, offvalue=0,
                                        height=1, width=8)

        ecpanel = tkinter.Frame(panel)
        eclabel = tkinter.Label(ecpanel, text="Endcond: ")
        endcondinput = ttk.Combobox(ecpanel, width=6,
                                    textvariable=self._endcondchoice)

        xinput["values"] = coords
        yinput["values"] = coords
        zinput["values"] = coords
        endcondinput["values"] = endconds

        xlabel.pack(side="left")
        xinput.pack(side="left")
        xltick.pack(side="left")
        ylabel.pack(side="left")
        yinput.pack(side="left")
        yltick.pack(side="left")
        zlabel.pack(side="left")
        zinput.pack(side="left")
        zltick.pack(side="left")

        eclabel.pack(side="left")
        endcondinput.pack(side="left")

        xpanel.pack()
        ypanel.pack()
        zpanel.pack()
        tickequal.pack()
        ecpanel.pack()

        self._plot()

    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        self._xcoord  = self._xchoice.get()
        self._ycoord  = self._ychoice.get()
        self._zcoord  = self._zchoice.get()
        self._equal   = self._equalchoice.get() == 1
        self._xlog    = self._xlogchoice.get()
        self._ylog    = self._ylogchoice.get()
        self._zlog    = self._zlogchoice.get()
        self._endcond = self._endcondchoice.get()

        if self._endcond == "all":
            self._endcond = None

        x = self._xcoord
        y = self._ycoord
        z = self._zcoord

        log = (self._xlog==1, self._ylog==1, self._zlog==1)

        if x == "None":
            x = None
        if y == "None":
            y = None

        if z == "None":
            ax = fig.add_subplot(1,1,1)
            self._orbits.plot(x, y, equal=self._equal, log=log,
                              endcond=self._endcond, axes=ax)
        else:
            # This would work for newer versions of matplotlib
            #ax = fig.add_subplot(1,1,1, projection="3d")
            ax = Axes3D(fig)
            self._orbits.plot(x, y, z, equal=self._equal, log=log,
                              endcond=self._endcond, axes=ax)

        self.draw()
