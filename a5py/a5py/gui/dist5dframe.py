"""
Contains definition of Dist5DFrame class.

File: dist5dframe.py
"""
import tkinter
import copy
import tkinter.ttk as ttk

from .plotframe import PlotFrame


class Dist5DFrame(PlotFrame):
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
                                      variable=self._plottype, value="1dexi",
                                      command=self._update_sidepanel)
        buttond = tkinter.Radiobutton(top, text="2D (E, xi)",
                                      variable=self._plottype, value="2dexi",
                                      command=self._update_sidepanel)

        buttona.pack(side="left")
        buttonb.pack(side="left")
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


    def _show_1dexipanel(self):
        """
        Show panel for plotting 1D E,xi dist.
        """
        data = self._dist.get_E_xi_dist()

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


    def _show_2dpanel(self):
        """
        Show panel for plotting 2D vpa,vpe dist.
        """
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
                                        height=1, width=8)

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

    def _show_2dexipanel(self):
        """
        Show panel for plotting 2D E,xi dist.
        """
        data = self._dist.get_E_xi_dist()

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
                                        height=1, width=8)

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
        """
        Set x coordinate for plots.

        Plot button is disabled if x and y coordinates are same.
        """
        self._xcoord = self._xchoice.get()

        if (self._plottype.get() == "2d" or self._plottype.get() == "2dexi") \
           and (self._xcoord == self._ycoord):
            self._plotbutton.config(state="disable")
        else:
            self._plotbutton.config(state="normal")

    def _set_ycoord(self, *args):
        """
        Set y coordinate for plots.

        Plot button is disabled if x and y coordinates are same.
        """
        self._ycoord = self._ychoice.get()

        if (self._plottype.get() == "2d" or self._plottype.get() == "2dexi") \
           and (self._xcoord == self._ycoord):
            self._plotbutton.config(state="disable")
        else:
            self._plotbutton.config(state="normal")

    def _plot(self):
        """
        Plot distribution whose view is displayed.
        """
        fig = self.get_fig()
        ax = fig.add_subplot(1,1,1)

        if(self._plottype.get() == "1d"):
            self._dist.plot_dist(self._xcoord, axes=ax)

        if(self._plottype.get() == "1dexi"):
            self._dist.plot_E_xi_dist(self._xcoord, axes=ax)

        elif(self._plottype.get() == "2d"):
            equal = self._equalaxis.get() == 1
            self._dist.plot_dist(self._xcoord, self._ycoord, equal=equal,
                                 axes=ax)

        elif(self._plottype.get() == "2dexi"):
            equal = self._equalaxis.get() == 1
            self._dist.plot_E_xi_dist(self._xcoord, self._ycoord, equal=equal,
                                      axes=ax)

        self.draw()
