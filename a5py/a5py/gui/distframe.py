"""
Contains definition of Dist5DFrame class.

File: dist5dframe.py
"""
import tkinter
import copy
import tkinter.ttk as ttk

from .plotframe import PlotFrame
from .components import NumEntry, DropdownMenu, Tickbox

class DistFrame(PlotFrame):
    """
    A frame for plotting 5D distribution interactively.

    The frame has several views which are toggled by radiobuttons on top panel.
    Each time a view is changed, all widgets in sidepanel are destroyed and new
    ones will be created.
    """

    def __init__(self, gui, dist):
        """
        Initialize frame and show default plot.

        Args:
            gui : GUI
            dist : Dist5D
        """
        super().__init__(gui)
        self._dist = dist

        self._xchoice   = tkinter.StringVar(self)
        self._ychoice   = tkinter.StringVar(self)
        self._xnum      = tkinter.IntVar(self)
        self._ynum      = tkinter.IntVar(self)
        self._equalaxis = tkinter.IntVar(self)
        self._logscale  = tkinter.IntVar(self)

        #self._xchoice.trace('w', self._change_coord)
        #self._ychoice.trace('w', self._change_coord)

        top = self.get_toppanel()

        self._plottype = tkinter.StringVar()
        buttona = tkinter.Radiobutton(top, text="1D (vpa, vpe)",
                                      variable=self._plottype, value="1d",
                                      command=self._update_sidepanel)
        buttonb = tkinter.Radiobutton(top, text="2D (vpa, vpe)",
                                      variable=self._plottype, value="2d",
                                      command=self._update_sidepanel)

        buttona.pack(side="left")
        buttonb.pack(side="left")

        # 6D dists have no E-xi dist yet.
        if hasattr(dist, "get_E_xi_dist"):
            buttonc = tkinter.Radiobutton(top, text="1D (E, xi)",
                                          variable=self._plottype,
                                          value="1dexi",
                                          command=self._update_sidepanel)
            buttond = tkinter.Radiobutton(top, text="2D (E, xi)",
                                          variable=self._plottype,
                                          value="2dexi",
                                          command=self._update_sidepanel)

            buttonc.pack(side="left")
            buttond.pack(side="left")

        self._plottype.set("2d")
        self._update_sidepanel()


    def _update_sidepanel(self):
        """
        Update side panel. Use this to change the view.
        """
        val = self._plottype.get()

        if val == "1d":
            self._show_1dpanel()
        elif val == "2d":
            self._show_2dpanel()
        elif val == "1dexi":
            self._show_1dexipanel()
        elif val == "2dexi":
            self._show_2dexipanel()

        self._plot()


    def _show_1dpanel(self):
        """
        Show panel for plotting 1D vpa,vpe dist.
        """
        data = self._dist.get_dist()

        coords = data["abscissae"]
        self._xcoord = coords[0]
        self._ycoord = None

        self._coords = copy.copy(coords)

        panel = self.get_sidepanel()

        self.plotbutton = tkinter.Button(panel, text="Plot",
                                         command=self._plot)

        self.xinput = DropdownMenu(panel, coords[0], self._coords, log=False,
                                   trace=None, label="x: ")
        self.ticklog = Tickbox(panel, 0, label="Log scale",
                               trace=self._plot)

        self.xinput.pack()
        self.ticklog.pack()
        self.plotbutton.pack()

        del data

        self._plot()


    def _show_1dexipanel(self):
        """
        Show panel for plotting 1D E,xi dist.
        """
        data = self._dist.get_E_xi_dist()

        coords = data["abscissae"]
        self._xcoord = coords[3]
        self._ycoord = None

        self._coords = copy.copy(coords)

        panel = self.get_sidepanel()

        self.plotbutton = tkinter.Button(panel, text="Plot",
                                         command=self._plot)

        self.xinput = DropdownMenu(panel, coords[3], self._coords, log=False,
                                   trace=None, label="x: ")
        self.Enum_entry  = NumEntry(panel, labeltext="E bins:",
                                    defval=10, isint=True)
        self.xinum_entry = NumEntry(panel, labeltext="xi bins:",
                                    defval=10, isint=True)

        self.ticklog = Tickbox(panel, 0, label="Log scale",
                               trace=self._plot)

        self.xinput.pack()
        self.Enum_entry.pack()
        self.xinum_entry.pack()
        self.ticklog.pack()
        self.plotbutton.pack()

        del data

        self._plot()


    def _show_2dpanel(self):
        """
        Show panel for plotting 2D vpa,vpe dist.
        """
        data = self._dist.get_dist()

        coords = data["abscissae"]
        self._xcoord = coords[0]
        self._ycoord = coords[2]

        self._coords = copy.copy(coords)

        panel = self.get_sidepanel()

        self.plotbutton = tkinter.Button(panel, text="Plot", command=self._plot)

        panel1 = tkinter.Frame(panel)
        self.xinput = DropdownMenu(panel1, coords[0], self._coords, log=False,
                                   trace=self._check_coord, label="x: ")
        self.yinput = DropdownMenu(panel1, coords[2], self._coords, log=False,
                                   trace=self._check_coord, label="y: ")

        panel2 = tkinter.Frame(panel)
        self.tickequal = Tickbox(panel2, 1, label="Axis equal",
                                 trace=self._plot)
        self.ticklog   = Tickbox(panel2, 0, label="Log scale",
                                 trace=self._plot)

        self.xinput.pack()
        self.yinput.pack()
        panel1.pack()

        self.tickequal.pack()
        self.ticklog.pack()
        panel2.pack()
        self.plotbutton.pack()

        del data

        self._plot()


    def _show_2dexipanel(self):
        """
        Show panel for plotting 2D E,xi dist.
        """
        data = self._dist.get_E_xi_dist()

        coords = data["abscissae"]
        self._xcoord = coords[0]
        self._ycoord = coords[2]

        self._coords = copy.copy(coords)

        panel = self.get_sidepanel()

        self.plotbutton = tkinter.Button(panel, text="Plot", command=self._plot)

        panel1 = tkinter.Frame(panel)
        self.xinput = DropdownMenu(panel1, coords[0], self._coords, log=False,
                                   trace=self._check_coord, label="x: ")
        self.yinput = DropdownMenu(panel1, coords[2], self._coords, log=False,
                                   trace=self._check_coord, label="y: ")

        panel2 = tkinter.Frame(panel)
        self.Enum_entry  = NumEntry(panel2, labeltext="E bins:",
                                    defval=10, isint=True)
        self.xinum_entry = NumEntry(panel2, labeltext="xi bins:",
                                    defval=10, isint=True)

        self.tickequal = Tickbox(panel2, 1, label="Axis equal",
                                 trace=self._plot)
        self.ticklog   = Tickbox(panel2, 0, label="Log scale",
                                 trace=self._plot)

        self.xinput.pack()
        self.yinput.pack()
        panel1.pack()

        self.Enum_entry.pack()
        self.xinum_entry.pack()
        self.tickequal.pack()
        self.ticklog.pack()
        panel2.pack()
        self.plotbutton.pack()

        del data

        self._plot()


    def _check_coord(self, *args):
        """
        Check that coordinates are valid.

        Plot button is disabled if x and y coordinates are same.
        """
        if (self._plottype.get() == "2d" or self._plottype.get() == "2dexi") \
           and (self.xinput.getval() == self.yinput.getval()):
            self.plotbutton.config(state="disable")
        else:
            self.plotbutton.config(state="normal")


    def _plot(self, *args):
        """
        Plot distribution whose view is displayed.
        """
        if (self._plottype.get() == "2d" or self._plottype.get() == "2dexi") \
           and (self.xinput.getval() == self.yinput.getval()):
            return

        fig = self.get_fig()
        ax = fig.add_subplot(1,1,1)
        logscale = self.ticklog.getval() == 1

        if(self._plottype.get() == "1d"):
            self._dist.plot_dist(self.xinput.getval(), logscale=logscale,
                                 axes=ax)

        if(self._plottype.get() == "1dexi"):
            E_edges  = self.Enum_entry.getval()
            xi_edges = self.xinum_entry.getval()
            self._dist.plot_E_xi_dist(self.xinput.getval(), logscale=logscale,
                                      E_edges=E_edges, xi_edges=xi_edges,
                                      axes=ax)

        elif(self._plottype.get() == "2d"):
            equal = self.tickequal.getval() == 1
            self._dist.plot_dist(self.xinput.getval(), self.yinput.getval(),
                                 equal=equal, logscale=logscale, axes=ax)

        elif(self._plottype.get() == "2dexi"):
            E_edges  = self.Enum_entry.getval()
            xi_edges = self.xinum_entry.getval()
            equal = self.tickequal.getval() == 1
            self._dist.plot_E_xi_dist(self.xinput.getval(),
                                      self.yinput.getval(),
                                      E_edges=E_edges, xi_edges=xi_edges,
                                      equal=equal, logscale=logscale, axes=ax)

        self.draw()
