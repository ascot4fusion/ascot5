"""
Contains definition of StateFrame class.

File: stateframe.py
"""
import copy
import tkinter
import tkinter.ttk as ttk

from mpl_toolkits.mplot3d import Axes3D

from .plotframe import PlotFrame
from .components import NumEntry, DropdownMenu, Tickbox
from .cameraControlPanel import generateCameraControlPanel



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
        self.avtk = None # So that we always regenerate the wall when entering here.
        if endstate is not None:
            self._states.append("endstate")

        # List of all possible coordinates.
        self._coords = ["R", "phimod", "z", "time", "energy", "pitch", "vnorm",
                        "Bnorm", "vR", "vphi", "vz", "BR", "Bphi", "Bz", "mu",
                        "vpar", "charge", "id", "x", "y", "phi", "rho",
                        "polmod", "pol", "None"]

        # All end conditions.
        if endstate is not None:
            endconds, counts = endstate.listendconds()
        else:
            endconds, counts = inistate.listendconds()
        self._endconds = ["all"] + endconds

        # Make top panel to choose what kind of plot is shown.
        top = self.get_toppanel()
        self.plotkind = tkinter.StringVar(self)
        self.plotkind.set("scatter")
        buttona = tkinter.Radiobutton(top, text="Scatter",
                                      value="scatter",
                                      variable=self.plotkind,
                                      command=self._showscatterpanel)
        buttonb = tkinter.Radiobutton(top, text="Histogram",
                                      value="histogram",
                                      variable=self.plotkind,
                                      command=self._showhistogrampanel)
        buttonc = tkinter.Radiobutton(top, text="3D wall-loads with VTK",
                                      value="3DwallLoad",
                                      variable=self.plotkind,
                                      command=self._show3dwallloadpanel)

        # Ini/endstate toggle
        self._statechoice = tkinter.StringVar(self)
        self._statechoice.set("inistate")
        if endstate is not  None:
            inibutton = tkinter.Radiobutton(top, text="Inistate",
                                            variable=self._statechoice,
                                            value="inistate",
                                            command=self._plot)
            endbutton = tkinter.Radiobutton(top, text="Endstate",
                                            variable=self._statechoice,
                                            value="endstate",
                                            command=self._plot)
            inibutton.grid(row=0, column=1)
            endbutton.grid(row=1, column=1)


        buttona.grid(row=0, column=0)
        buttonb.grid(row=1, column=0)
        buttonc.grid(row=2, column=0)

        self._showscatterpanel()

    def _show3dwallloadpanel(self):
        self._statechoice.set("endstate")
        panel = self.get_sidepanel()
        __fig = self.get_fig() #clears the plots, at least...
        
        #self.caxpanel   = tkinter.Frame(panel)
        self.Pmin_entry = NumEntry(panel, labeltext="P [W/m2] min:", defval=100.0)
        self.Pmax_entry = NumEntry(panel, labeltext="P [W/m2] max:", defval=1.0e6)
        self.Plogtick   = Tickbox(panel, 1, label="log10 P (color)")

        #self.Pmin_entry.grid(row=0, column=0, sticky="W")
        #self.Pmax_entry.grid(row=1, column=0, sticky="W")
        #self.Plogtick.grid(  row=2, column=0, sticky="W")
        self.Pmin_entry.pack()
        self.Pmax_entry.pack()
        self.Plogtick.pack()
        
        self.camControlPanelComponents,self.applyCameraControlPanel= generateCameraControlPanel(panel)
        self.camControlPanelComponents['panel'].pack()


        self.PplotBtn   = tkinter.Button(panel, text="Plot with VTK", command=self._plotVTK)
        
        #self.caxpanel.pack()
        #self.PplotBtn.grid(  row=3, column=0, sticky="W")
        self.PplotBtn.pack()

        
        self.draw()

    def _plotVTK(self):
        print('Importing VTK.')
        import a5py.wallloads.toVTK
        import a5py.wall.a5vtkwall
        import numpy as np

        print('Starting to plot with VTK.')
        
        if self.avtk is None:
            ASCOT = self._gui.get_ascotobject()
            self.avtk = a5py.wallloads.toVTK.as3DpolyData(ASCOT)

        colorRange=(np.log10(float(self.Pmin_entry.getval())),
                    np.log10(float(self.Pmax_entry.getval()))  )
        
        logScale = self.Plogtick.getval() == 1
        
        self.camControl = a5py.wall.a5vtkwall.camControl()
        
        self.applyCameraControlPanel(self.camControlPanelComponents,self.camControl)
        
        self.avtk.plot(manual_range=colorRange,
                       logarithmicColorScale=logScale,
                       camControl=self.camControl)


    def _showscatterpanel(self):
        panel = self.get_sidepanel()

        self.xchoice = DropdownMenu(panel, "R", self._coords, log=True,
                                    trace=self._plotscatter, label="x: ")
        self.ychoice = DropdownMenu(panel, "z", self._coords, log=True,
                                    trace=self._plotscatter, label="y: ")
        self.zchoice = DropdownMenu(panel, "None", self._coords, log=True,
                                    trace=self._plotscatter, label="z: ")
        self.cchoice = DropdownMenu(panel, "None", self._coords, log=True,
                                    trace=self._plotscatter, label="c: ")

        self.endcondchoice = DropdownMenu(panel, "all", self._endconds,
                                          trace=self._plotscatter,
                                          label="Endcond: ")

        self.equalchoice = Tickbox(panel, 1, label="Axis equal",
                                   trace=self._plotscatter)

        self.xchoice.pack()
        self.ychoice.pack()
        self.zchoice.pack()
        self.cchoice.pack()
        self.equalchoice.pack()
        self.endcondchoice.pack()

        self._plotscatter()


    def _showhistogrampanel(self):
        panel = self.get_sidepanel()

        self.xchoice = DropdownMenu(panel, "R", self._coords, log=True,
                                    trace=self._plothistogram, label="x: ")
        self.ychoice = DropdownMenu(panel, "None", self._coords, log=True,
                                    trace=self._plothistogram, label="y: ")

        self.zlogtick = Tickbox(panel, 1, label="log10 z",
                                trace=self._plothistogram)
        self.weighttick = Tickbox(panel, 1, label="weight",
                                  trace=self._plothistogram)

        self.endcondchoice = DropdownMenu(panel, "all", self._endconds,
                                          trace=self._plothistogram,
                                          label="Endcond: ")

        binpanel   = tkinter.Frame(panel)
        self.xmin_entry = NumEntry(binpanel, labeltext="x_min:", defval=0)
        self.xmax_entry = NumEntry(binpanel, labeltext="x_max:", defval=10)
        self.xnum_entry = NumEntry(binpanel, labeltext="x_num:", defval=10,
                                   isint=True)

        self.ymin_entry = NumEntry(binpanel, labeltext="y_min:", defval=1)
        self.ymax_entry = NumEntry(binpanel, labeltext="y_max:", defval=2)
        self.ynum_entry = NumEntry(binpanel, labeltext="y_num:", defval=10,
                                   isint=True)

        self.xmin_entry.grid(row=0, column=0, sticky="W")
        self.xmax_entry.grid(row=1, column=0, sticky="W")
        self.xnum_entry.grid(row=2, column=0, sticky="W")

        self.ymin_entry.grid(row=0, column=1, sticky="W")
        self.ymax_entry.grid(row=1, column=1, sticky="W")
        self.ynum_entry.grid(row=2, column=1, sticky="W")

        self.xchoice.pack()
        self.ychoice.pack()
        self.zlogtick.pack()
        self.weighttick.pack()
        self.endcondchoice.pack()
        binpanel.pack()
        tkinter.Button(panel, text="Plot", command=self._plothistogram).pack()

        self._plothistogram()


    def _plotscatter(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        if self._statechoice.get() == "inistate":
            state = self._inistate
        else:
            state = self._endstate

        xcoord  = self.xchoice.getval()
        ycoord  = self.ychoice.getval()
        zcoord  = self.zchoice.getval()
        ccoord  = self.cchoice.getval()
        equal   = self.equalchoice.getval() == 1
        xlog    = self.xchoice.islog()
        ylog    = self.ychoice.islog()
        zlog    = self.zchoice.islog()
        clog    = self.cchoice.islog()
        endcond = self.endcondchoice.getval()

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
            axes = fig.add_subplot(1,1,1, projection="3d")
            state.scatter(xcoord, ycoord, zcoord, ccoord,
                          equal=equal, log=log,
                          endcond=endcond, axes=axes)

        self.draw()


    def _plothistogram(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        if self._statechoice.get() == "inistate":
            state = self._inistate
        else:
            state = self._endstate

        xcoord  = self.xchoice.getval()
        ycoord  = self.ychoice.getval()
        if xcoord == "None":
            xcoord = None
        if ycoord == "None":
            ycoord = None

        endcond = self.endcondchoice.getval()
        if endcond == "all":
            endcond = None

        logx     = self.xchoice.islog()
        logy     = self.ychoice.islog()
        logscale = self.zlogtick.getval() == 1

        xbins = [ self.xmin_entry.getval(),
                  self.xmax_entry.getval(),
                  self.xnum_entry.getval() + 1 ]

        ybins = [ self.ymin_entry.getval(),
                  self.ymax_entry.getval(),
                  self.ynum_entry.getval() + 1 ]

        weight = self.weighttick.getval() == 1

        axes = fig.add_subplot(1,1,1)
        state.histogram(x=xcoord, y=ycoord, xbins=xbins, ybins=ybins,
                        weight=weight, logx=logx, logy=logy,
                        logscale=logscale, endcond=endcond, axes=axes)

        self.draw()


    def _plot(self):
        if self.plotkind.get() == "scatter":
            self._plotscatter()
        else:
            self._plothistogram()


