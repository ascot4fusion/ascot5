"""
State plotting (frame) for the output content.

File: outputstate.py
"""
import tkinter as tk
from tkinter import ttk
import numpy as np

from collections import OrderedDict
from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox, ToggleButton

class StateFrame(ttk.Frame):
    """
    Settings frame for visualizing ini/end state.

    Visualizing is done either via scatter plot (2D/3D + color) or histogram
    (1D/2D).

    This frame contains Notebook widget containing two nested frames: one for
    the scatter plot settings and the other for the histogram settings. Both
    share the same StateCanvas frame.
    """

    def init(self, gui, canvas):
        """
        Define nested frames and initialize the notebook.
        """
        self.gui = gui

        # Quantities that can be plotted (should be generated from state.py)
        quantities = OrderedDict([
            ("R",         "m"),
            ("phi (mod)", "deg"),
            ("z",         "m"),
            ("Time",      "s"),
            ("Energy",    "eV"),
            ("Pitch",     "vpa/|v|"),
            ("ID",        "1"),
            ("x",         "m"),
            ("y",         "m"),
            ("rho",       "1"),
            ("pol (mod)", "deg"),
            ("phi",       "deg"),
            ("pol",       "deg"),
            ("mu",        "eV/T"),
            ("|p|",       "kg m/s"),
            ("ppa",       "kg m/s"),
            ("ppe",       "kg m/s"),
            ("pR",        "kg m/s"),
            ("pphi",      "kg m/s"),
            ("pz",        "kg m/s"),
            ("|B|",       "T"),
            ("BR",        "T"),
            ("Bphi",      "T"),
            ("Bz",        "T"),
            ("charge",    "T")
        ])

        # Define separate classes for the nested frames
        class ScatterFrame(ttk.Frame):
            """
            Settings for the scatter plot.
            """

            def init(self, canvas):
                """
                Initialize all the widgets.
                """

                # For the overall layout we use three frames
                f1 = ttk.Frame(self) # For coordinate selection widgets
                f2 = ttk.Frame(self) # Plot and store buttons
                f3 = ttk.Frame(self) # Additional settings

                f1.grid(row=0, column=0, sticky="e")
                f2.grid(row=0, column=1, sticky="ne")
                f3.grid(row=1, column=0, sticky="ew")

                self.columnconfigure(1, weight=1)

                # Choose coordinate for [x,y,z,c] and tick if log scale is used
                self.xcrd = DropdownMenu(f1, log=True, label="x: ", width=12)
                self.ycrd = DropdownMenu(f1, log=True, label="y: ", width=12)
                self.zcrd = DropdownMenu(f1, log=True, label="z: ", width=12)
                self.ccrd = DropdownMenu(f1, log=True, label="c: ", width=12)

                self.xcrd.grid(row=0, column=0)
                self.ycrd.grid(row=1, column=0)
                self.zcrd.grid(row=2, column=0)
                self.ccrd.grid(row=3, column=0)

                # Is coordinate from ini or endstate
                self.xbtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.ybtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.zbtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.cbtn = ToggleButton(f1, label1text="Ini", label2text="End")

                self.xbtn.grid(row=0, column=1)
                self.ybtn.grid(row=1, column=1)
                self.zbtn.grid(row=2, column=1)
                self.cbtn.grid(row=3, column=1)

                # Plot and store buttons
                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                # Only markers with this endstate are plotted
                self.endc = DropdownMenu(f3, label="Endcond: ", width=12,
                                         labelwidth=8)
                # Set x,y,z axes equal aspect ration
                self.axeq = Tickbox(f3, label=" Axis equal")

                self.endc.pack(side="left", anchor="w")
                self.axeq.pack(side="left", anchor="w")

                # Set values for the coordinate drop-down menus
                self.xcrd.setvals(list(quantities.keys()), "R", 0)
                self.ycrd.setvals(list(quantities.keys()), "z", 0)
                # For z and c one can choose "None" which is also the default
                qwithnone = quantities.copy()
                qwithnone["None"] = ("None","")
                self.zcrd.setvals(list(qwithnone.keys()), "None", 0)
                self.ccrd.setvals(list(qwithnone.keys()), "None", 0)

                self.canvas = canvas
                return self

            def plot(self, run):
                """
                Plot the given run using the active settings in the frame.
                """

                # Read and process settings
                xcoord  = self.xcrd.getval()
                ycoord  = self.ycrd.getval()
                zcoord  = self.zcrd.getval()
                ccoord  = self.ccrd.getval()
                equal   = self.axeq.getval() == 1
                xlog    = self.xcrd.islog()
                ylog    = self.ycrd.islog()
                zlog    = self.zcrd.islog()
                clog    = self.ccrd.islog()
                endcond = self.endc.getval()

                if endcond == "all":
                    endcond = None

                log = (xlog==1, ylog==1, zlog==1, clog==1)

                if ccoord == "None":
                    ccoord = None

                iniend = ["e" if self.xbtn.var.get() else "i",
                          "e" if self.ybtn.var.get() else "i",
                          "e" if self.zbtn.var.get() else "i",
                          "e" if self.cbtn.var.get() else "i"]

                if zcoord == "None":
                    # 2D plot
                    self.canvas.clear()
                    run.plotstate_scatter(xcoord, ycoord, c=ccoord,
                                          axesequal=equal, log=log,
                                          endcond=endcond, iniend=iniend,
                                          axes=self.canvas.axes)
                else:
                    # 3D plot (here we have to first create 3D axes for canvas)
                    self.canvas.make3d()
                    self.canvas.clear()
                    run.plotstate_scatter(xcoord, ycoord, zcoord, ccoord,
                                          axesequal=equal, log=log,
                                          iniend=iniend, endcond=endcond,
                                          axes=self.canvas.axes)
                self.canvas.draw()


        class HistFrame(ttk.Frame):
            """
            Settings for the histogram frame.
            """

            def init(self, canvas):
                """
                Initialize widgets.
                """

                # Use four frames for layout
                f1 = ttk.Frame(self) # Coordinate selection
                f2 = ttk.Frame(self) # Plot and save buttons
                f3 = ttk.Frame(self) # Bin selection
                f4 = ttk.Frame(self) # Other settings

                f1.grid(row=0, column=0, sticky="e")
                f2.grid(row=0, column=1, sticky="ne")
                f3.grid(row=1, column=0, columnspan=2, sticky="w")
                f4.grid(row=2, column=0, sticky="w")

                self.columnconfigure(1, weight=1)

                # Coordinate selection
                self.xcrd = DropdownMenu(f1, log=True, label="x: ", width=18,
                                         trace=self.newxcoord)
                self.ycrd = DropdownMenu(f1, log=True, label="y: ", width=18,
                                         trace=self.newycoord)

                # Coordinate taken from ini or endstate
                self.xbtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.ybtn = ToggleButton(f1, label1text="Ini", label2text="End")

                self.xcrd.grid(row=0, column=0)
                self.ycrd.grid(row=1, column=0)

                self.xbtn.grid(row=0, column=1)
                self.ybtn.grid(row=1, column=1)

                # Plot and save buttons
                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                # Weighted histogram (or no), axes equal aspect ratio (or no)
                self.wght = Tickbox(f4, width=1)
                self.axeq = Tickbox(f4, width=1)
                self.zlog = Tickbox(f4, width=1)

                self.wght.grid(row=0,column=0,sticky="w")
                tk.Label(f4, text="With weights").grid(row=0,column=1,sticky="w")
                self.axeq.grid(row=1,column=0,sticky="w")
                tk.Label(f4, text="Axis equal").grid(row=1,column=1,sticky="w")
                self.zlog.grid(row=2,column=0,sticky="w")
                tk.Label(f4, text="log10 scale").grid(row=2,column=1,sticky="w")

                # Set NumEntries for xmin, xmax, and nx inputs
                self.xmin_entry = NumEntry(
                    f3, labeltext="rho [1] = ", entrywidth=5, labelwidth=22,
                    anchor="e", defval=0.0)
                self.xmax_entry = NumEntry(f3, labeltext="–", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=1.0)
                self.xnum_entry = NumEntry(f3, labeltext="x", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=50,
                                           isint=True)

                # Set NumEntries for ymin, ymax, and ny inputs
                self.ymin_entry = NumEntry(
                    f3, labeltext="* = ", entrywidth=5, labelwidth=22,
                    anchor="e", defval=0.0)
                self.ymax_entry = NumEntry(f3, labeltext="–", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=1.0)
                self.ynum_entry = NumEntry(f3, labeltext="x", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=50,
                                           isint=True)


                self.xmin_entry.grid(row=0, column=0)
                self.xmax_entry.grid(row=0, column=1)
                self.xnum_entry.grid(row=0, column=2)
                self.ymin_entry.grid(row=1, column=0)
                self.ymax_entry.grid(row=1, column=1)
                self.ynum_entry.grid(row=1, column=2)

                # Set values for the drop-down menus and allow y coordinate to
                # have None option (in which case the histogram is 1D)
                self.xcrd.setvals(list(quantities.keys()), "rho", 0)
                qwithnone = quantities.copy()
                qwithnone["None"] = ("None","")
                self.ycrd.setvals(list(qwithnone.keys()), "None", 0)

                #gui.params.add(
                    #output_stateplot_minr=xmin_entry.choice,
                    #output_stateplot_maxr=xmax_entry.choice,
                    #output_stateplot_numr=xnum_entry.choice,
                    #output_stateplot_minz=ymin_entry.choice,
                    #output_stateplot_maxz=ymax_entry.choice,
                    #output_stateplot_numz=ynum_entry.choice,
                    #output_stateplot_qnt=xcrd.var
                    #)
                self.canvas = canvas

                return self


            def newxcoord(self, *args):
                """
                x-coordinate changed, update the label on the x bin selection.
                """
                val = self.xcrd.var.get()
                unit = quantities[val]
                self.xmin_entry.setlabel(val + " [" + unit + "]" + " = ")


            def newycoord(self, *args):
                """
                y-coordinate changed, update the label on the y bin selection.

                If None, we disable the bin selection completely for y.
                """
                val = self.ycrd.var.get()
                if val == "None":
                    self.ymin_entry.disable()
                    self.ymax_entry.disable()
                    self.ynum_entry.disable()
                    self.ymin_entry.setlabel("* = ")
                else:
                    unit = quantities[val]
                    self.ymin_entry.enable()
                    self.ymax_entry.enable()
                    self.ynum_entry.enable()
                    self.ymin_entry.setlabel(val + " [" + unit + "]" + " = ")


            def plot(self, run):
                """
                Plot the histogram for the given run using active settings.
                """
                xcoord = self.xcrd.getval()
                ycoord = self.ycrd.getval()
                if ycoord == "None": ycoord = None

                logx      = self.xcrd.islog()
                logy      = self.ycrd.islog()
                logscale  = self.zlog.getval() == 1
                axesequal = self.axeq.getval() == 1
                weight    = self.wght.getval() == 1

                iniend = ["e" if self.xbtn.var.get() else "i",
                          "e" if self.ybtn.var.get() else "i"]

                def getbins(vmin, vmax, nv, islog):
                    """
                    Short function to check user input and generate bins.
                    """
                    if vmax < vmin or nv <= 0: return None
                    if islog and (vmin <= 0 or vmax <= 0): return None

                    if islog:
                        return np.logspace(np.log10(vmin), np.log10(vmax), nv)
                    else:
                        return np.linspace(vmin, vmax, nv)

                xbins = getbins(self.xmin_entry.getval(),
                                self.xmax_entry.getval(),
                                self.xnum_entry.getval(),
                                logx)
                ybins = getbins(self.ymin_entry.getval(),
                                self.ymax_entry.getval(),
                                self.ynum_entry.getval(),
                                logy)

                self.canvas.clear()
                run.plotstate_histogram(
                    xcoord, y=ycoord, xbins=xbins, ybins=ybins,
                    endcond=None, weight=weight, iniend=iniend,
                    log=[logx, logy], logscale=logscale,
                    axesequal=axesequal, axes=self.canvas.axes)
                self.canvas.draw()

        # Nested frames defined. Now we just create instances and put the onto
        # the notebook located in this main frame.
        master = ttk.Notebook(self)

        self.framescatter = ScatterFrame(master).init(canvas)
        self.framehist    = HistFrame(master).init(canvas)
        master.add(self.framescatter, text="Scatter")
        master.add(self.framehist, text="Histogram")
        master.pack(fill="both", expand=True)

        self.canvas = canvas
        return self


    def display(self):
        """
        Display state frame for the run that is currently active.
        """
        run = self.gui.ascot.hdf5.active
        try:
            endconds, counts = run.endstate.listendconds()
            if self.framescatter.xbtn.switch.isdisabled:
                self.framescatter.xbtn.enable()
                self.framescatter.ybtn.enable()
                self.framescatter.zbtn.enable()
                self.framescatter.cbtn.enable()
                self.framehist.xbtn.enable()
                self.framehist.ybtn.enable()

        except AttributeError:
            endconds, counts = run.inistate.listendconds()
            self.framescatter.xbtn.disable()
            self.framescatter.ybtn.disable()
            self.framescatter.zbtn.disable()
            self.framescatter.cbtn.disable()
            self.framehist.xbtn.disable()
            self.framehist.ybtn.disable()

        endconds = ["all"] + endconds

        self.framescatter.endc.setvals(endconds, "all")

        def plotscatter():
            self.framescatter.plot(run)

        def plothist():
            self.framehist.plot(run)

        self.framescatter.plotbutton.configure(command=plotscatter)
        self.framehist.plotbutton.configure(command=plothist)


class StateCanvas(PlotFrame):
    """
    Canvas for plotting the state plots.
    """

    def __init__(self, master):
        """
        Add attribute axes3d which flags whether next axes should be 3d.
        """
        self.axes3d = False
        super().__init__(master)

    def make3d(self):
        """
        Flag next axes to be created 3d.
        """
        self.axes3d = True

    def set_axes(self):
        """
        Override this method to allow generation of 3d axes.
        """
        if self.axes3d:
            self.axes3d = False
            return self.fig.add_subplot(1,1,1, projection="3d")

        return self.fig.add_subplot(1,1,1)

    def init(self):
        """
        Do nothing.
        """
        return self
