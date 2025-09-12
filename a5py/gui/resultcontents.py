import tkinter as tk
from tkinter import ttk
import numpy as np

from matplotlib.gridspec import GridSpec
from .components import ContentTab, PlotFrame, NumEntry, DropdownMenu, Tickbox,\
    ToggleButton
from a5py.plotting.plotting import defaultcamera

class Summary(ContentTab):
    """Settings frame summarizing simulation output.

    The summary canvas is not interactive, so this frame does not contain any
    settings. Instead, it has a text box that shows sum of markers by
    end condition and possible error messages.
    """

    def __init__(self, frame, canvas, gui):
        """Initializes frame by creating all the widgets.
        """
        super().__init__(frame, canvas, gui)

        self.plot = Summary.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)
        self.markersummary = tk.Text(self, width=50)
        self.markersummary.pack(fill="both", expand=True)
        self.markersummary.config(state="disabled")

    def selecttab(self):
        """Display this frame and redraws canvas.
        """
        run = self.gui.ascot.data.active
        self.markersummary.config(state="normal")
        self.markersummary.delete("1.0", "end")

        text = "Marker end conditions:\n"
        end, err = run.getstate_markersummary()
        for m in end:
            text += m[1] + " : " + str(m[0]) + "\n"
        text += "\nErrors:\n"
        for m in err:
            text += m[0] + " at line " + str(m[1]) + " in file " + m[2] + "\n"

        self.markersummary.insert("end", text)
        self.markersummary.config(state="disabled")
        self.plot.clear()

        self.plot.axes_inirho.set_xlim(0,1.1)
        self.plot.axes_inirho.set_title("Initial radial position")
        run.plotstate_histogram("ini rho", xbins=np.linspace(0,1.1,55),
                                weight=True, axes=self.plot.axes_inirho)

        self.plot.axes_endrho.set_xlim([0,1.1])
        self.plot.axes_endrho.set_title("Final radial position")
        run.plotstate_histogram("end rho", xbins=np.linspace(0,1.1,55),
                                weight=True, axes=self.plot.axes_endrho)

        self.plot.axes_mileage.set_title("Final mileage")
        run.plotstate_histogram("log end mileage", xbins=55, weight=True,
                                axes=self.plot.axes_mileage)

        self.plot.axes_energy.set_title("Final energy")
        run.plotstate_histogram("log ini ekin", xbins=55, weight=True,
                                axes=self.plot.axes_energy)

        self.plot.axes_rz.set_title("Final R-z positions")
        run.plotstate_scatter("end R", "end z", c="C0", endcond=None,
                              axesequal=True, axes=self.plot.axes_rz)
        #self.gui.ascot.input_plotseparatrix(0, 0, axes=self.plot.axes_rz)
        self.gui.ascot.input_plotwallcontour(axes=self.plot.axes_rz)

        self.plot.axes_rhophi.set_xlim([0,1.2])
        self.plot.axes_rhophi.set_ylim([0,360])
        self.plot.axes_rhophi.set_xticks([0, 0.5, 1.0])
        self.plot.axes_rhophi.set_yticks([0, 180, 360])
        self.plot.axes_rhophi.set_title("Final rho-phi positions")
        run.plotstate_scatter("end rho", "end phimod", c="C0", endcond=None,
                              axesequal=False, axes=self.plot.axes_rhophi)

        self.plot.draw()
        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        """Canvas for showing summary plots.

        This canvas creates several axes where the plots can be placed.
        """

        def __init__(self, frame):
            """Name axes just so that they are easier to reference
            when plotting.
            """
            super().__init__(frame)
            self.axes_inirho  = None
            self.axes_endrho  = None
            self.axes_mileage = None
            self.axes_energy  = None
            self.axes_rz      = None
            self.axes_rhophi  = None

        def set_axes(self):
            """Override the PlotFrame method in order to create multiple
            axes.
            """
            gs = GridSpec(2,4)
            self.axes_inirho  = self.fig.add_subplot(gs[0,0])
            self.axes_endrho  = self.fig.add_subplot(gs[0,1])
            self.axes_mileage = self.fig.add_subplot(gs[1,0])
            self.axes_energy  = self.fig.add_subplot(gs[1,1])
            self.axes_rz      = self.fig.add_subplot(gs[:,2])
            self.axes_rhophi  = self.fig.add_subplot(gs[:,3])

            ax = [self.axes_inirho, self.axes_endrho, self.axes_mileage,
                  self.axes_energy, self.axes_rz, self.axes_rhophi]
            return ax

class StateScatter(ContentTab):
    """Settings frame for visualizing ini/end state.

    Visualizing is done either via scatter plot (2D/3D + color) or histogram
    (1D/2D).

    This frame contains Notebook widget containing two nested frames: one for
    the scatter plot settings and the other for the histogram settings. Both
    share the same StateCanvas frame.
    """

    def __init__(self, frame, canvas, gui):
        """Define nested frames and initialize the notebook.
        """
        super().__init__(frame, canvas, gui)
        self.plot = StateScatter.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

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

        def plot():
            """Plot the given run using the active settings in the frame.
            """
            # Read and process settings
            equal   = self.axeq.getval() == 1
            endcond = self.endc.getval()

            if endcond == "all":
                endcond = None

            coords = [None] * 4
            for i, coord in enumerate(["x", "y", "z", "c"]):
                coords[i] = getattr(self, coord + "crd").getval()
                if coords[i] == "None":
                    coords[i] = None
                    continue

                log = "log" if getattr(self, coord + "crd").islog() else ""
                ini = "end" if getattr(self, coord + "btn").var.get() else "ini"
                coords[i] = log + " " + ini + " " + coords[i]

            run = self.gui.ascot.data.active
            if coords[2] is None:
                # 2D plot
                self.plot.clear()
                run.plotstate_scatter(coords[0], coords[1], c=coords[3],
                                      axesequal=equal, endcond=endcond,
                                      axes=self.plot.axes)
            else:
                # 3D plot (here we have to first create 3D axes for canvas)
                self.plot.make3d()
                self.plot.clear()
                run.plotstate_scatter(coords[0], coords[1], coords[2],
                                      c=coords[3], axesequal=equal,
                                      endcond=endcond, axes=self.plot.axes)
            self.plot.draw()

        self.plotbutton.configure(command=plot)
        #self.savebutton.configure(command=save)

    def selecttab(self):
        run = self.gui.ascot.data.active

        # Set values for the coordinate drop-down menus
        quantities = run.getstate_list()
        self.xcrd.setvals(list(quantities.keys()), "r", 0)
        self.ycrd.setvals(list(quantities.keys()), "z", 0)

        # For z and c one can choose "None" which is also the default
        qwithnone = quantities.copy()
        qwithnone["None"] = ("None","")
        self.zcrd.setvals(list(qwithnone.keys()), "None", 0)
        self.ccrd.setvals(list(qwithnone.keys()), "None", 0)

        end, err = run.getstate_markersummary()
        _, end = zip(*end)
        end = ["all"] + list(end)
        self.endc.setvals(end, "all")

        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        """Canvas for plotting the state plots.
        """

        def __init__(self, canvas):
            """Add attribute axes3d which flags whether next axes should be 3d.
            """
            self.axes3d = False
            super().__init__(canvas)

        def make3d(self):
            """Flag next axes to be created 3d.
            """
            self.axes3d = True

        def set_axes(self):
            """Override this method to allow generation of 3d axes.
            """
            if self.axes3d:
                self.axes3d = False
                return self.fig.add_subplot(1,1,1, projection="3d")

            return self.fig.add_subplot(1,1,1)

class StateHistogram(ContentTab):
    """Settings for the histogram frame.
    """

    def __init__(self, frame, canvas, gui):
        """Initialize widgets.
        """
        super().__init__(frame, canvas, gui)
        self.plot = StateHistogram.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

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

        #gui.params.add(
        #output_stateplot_minr=xmin_entry.choice,
        #output_stateplot_maxr=xmax_entry.choice,
        #output_stateplot_numr=xnum_entry.choice,
        #output_stateplot_minz=ymin_entry.choice,
        #output_stateplot_maxz=ymax_entry.choice,
        #output_stateplot_numz=ynum_entry.choice,
        #output_stateplot_qnt=xcrd.var
        #)

        def plot():
            """Plot the histogram for the given run using active settings.
            """
            xcoord = self.xcrd.getval()
            xcoord = "end " + xcoord if self.xbtn.var.get() else "ini " + xcoord
            if self.xcrd.islog(): xcoord = "log " + xcoord

            ycoord = self.ycrd.getval()
            if ycoord == "None":
                ycoord = None
            else:
                ycoord = "end " + ycoord if self.ybtn.var.get() \
                    else "ini " + ycoord
                if self.ycrd.islog(): ycoord = "log " + ycoord

            logscale  = self.zlog.getval() == 1
            axesequal = self.axeq.getval() == 1
            weight    = self.wght.getval() == 1

            def getbins(vmin, vmax, nv, islog):
                """Short function to check user input and generate bins.
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
                            False)
            ybins = getbins(self.ymin_entry.getval(),
                            self.ymax_entry.getval(),
                            self.ynum_entry.getval(),
                            False)

            self.plot.clear()
            run = self.gui.ascot.data.active
            run.plotstate_histogram(
                xcoord, y=ycoord, xbins=xbins, ybins=ybins,
                endcond=None, weight=weight, logscale=logscale,
                axesequal=axesequal, axes=self.plot.axes)
            self.plot.draw()

        self.plotbutton.configure(command=plot)
        #self.savebutton.configure(command=save)

    def newxcoord(self, *args):
        """x-coordinate changed, update the label on the x bin selection.
        """
        val  = self.xcrd.var.get()
        try:
            unit = str(self.gui.ascot.data.active.getstate(val).units)
        except AttributeError:
            unit = "1"
        self.xmin_entry.setlabel(val + " [" + unit + "]" + " = ")


    def newycoord(self, *args):
        """y-coordinate changed, update the label on the y bin selection.

        If None, we disable the bin selection completely for y.
        """
        val = self.ycrd.var.get()
        if val == "None":
            self.ymin_entry.disable()
            self.ymax_entry.disable()
            self.ynum_entry.disable()
            self.ymin_entry.setlabel("* = ")
        else:
            try:
                unit = str(self.gui.ascot.data.active.getstate(val).units)
            except AttributeError:
                unit = "1"
            self.ymin_entry.enable()
            self.ymax_entry.enable()
            self.ynum_entry.enable()
            self.ymin_entry.setlabel(val + " [" + unit + "]" + " = ")

    def selecttab(self):
        """Display state frame for the run that is currently active.
        """
        run = self.gui.ascot.data.active
        # Set values for the drop-down menus and allow y coordinate to
        # have None option (in which case the histogram is 1D)
        quantities = run.getstate_list()
        self.xcrd.setvals(list(quantities.keys()), "rho", 0)
        qwithnone = quantities.copy()
        qwithnone["None"] = ("None","")
        self.ycrd.setvals(list(qwithnone.keys()), "None", 0)
        try:
            if self.xbtn.switch.isdisabled:
                self.xbtn.enable()
                self.ybtn.enable()
                self.zbtn.enable()
                self.cbtn.enable()
                self.xbtn.enable()
                self.ybtn.enable()

        except AttributeError:
            self.xbtn.disable()
            self.ybtn.disable()
            self.zbtn.disable()
            self.cbtn.disable()
            self.xbtn.disable()
            self.ybtn.disable()

        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        """Canvas for plotting the state plots.
        """
        pass

class Orbit(ContentTab):

    def __init__(self, frame, canvas, gui):
        """Initializes frame by creating all the widgets.
        """
        super().__init__(frame, canvas, gui)
        self.plot = Orbit.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)
        f3 = ttk.Frame(self)

        f1.grid(row=0, column=0, sticky="w")
        f2.grid(row=0, column=1, sticky="ne")
        f3.grid(row=1, column=0, sticky="ew")

        self.columnconfigure(1, weight=1)

        self.xcrd = DropdownMenu(f1, log=True, label="x: ", width=12,
                                 labelwidth=2)
        self.ycrd = DropdownMenu(f1, log=True, label="y: ", width=12,
                                 labelwidth=2)
        self.zcrd = DropdownMenu(f1, log=True, label="z: ", width=12,
                                 labelwidth=2)
        self.ccrd = DropdownMenu(f1, log=True, label="c: ", width=12,
                                 labelwidth=2)

        self.xcrd.grid(row=0, column=0)
        self.ycrd.grid(row=1, column=0)
        self.zcrd.grid(row=2, column=0)
        self.ccrd.grid(row=3, column=0)

        self.plotbutton = tk.Button(f2, text="Plot", width=3)
        self.savebutton = tk.Button(f2, text="Store", width=3)
        self.plotbutton.pack(anchor="nw")
        self.savebutton.pack(anchor="nw")

        self.endc = DropdownMenu(f3, label="Endcond: ", width=12,
                                 labelwidth=8)
        self.axeq = Tickbox(f3, label=" Axes equal")

        self.endc.pack(side="left", anchor="w")
        self.axeq.pack(side="left", anchor="w")

        def plot():
            equal   = self.axeq.getval() == 1
            endcond = self.endc.getval()

            if endcond == "all":
                endcond = None

            coords = [None] * 4
            for i, coord in enumerate(["x", "y", "z", "c"]):
                coords[i] = getattr(self, coord + "crd").getval()
                if coords[i] == "None":
                    coords[i] = None
                    continue

                if getattr(self, coord + "crd").islog():
                    coords[i] = "log " + coords[i]

            run = self.gui.ascot.data.active
            if coords[2] is None:
                self.plot.clear()
                run.plotorbit_trajectory(coords[0], coords[1], c=coords[3],
                                         axesequal=equal, endcond=endcond,
                                         axes=self.plot.axes)
                #if ccoord == "None":
                #    run.orbit.plot(xcoord, ycoord, equal=equal, log=log,
                #                   endcond=endcond, axes=self.plot.axes)
                #else:
                #    run.orbit.scatter(xcoord, ycoord, c=ccoord, equal=equal,
                #                      endcond=endcond, axes=self.plot.axes)
            else:
                self.plot.make3d()
                self.plot.clear()
                run.plotorbit_trajectory(coords[0], coords[1], coords[2],
                                         c=coords[3], axesequal=equal,
                                         endcond=endcond, axes=self.plot.axes)
                #if ccoord == "None":
                #    run.orbit.plot(xcoord, ycoord, zcoord, equal=equal,
                #                   endcond=endcond, axes=self.plot.axes)
                #else:
                #    run.orbit.scatter(xcoord, ycoord, zcoord, ccoord,
                #                      equal=equal, log=log,
                #                      endcond=endcond, axes=self.plot.axes)

            self.plot.draw()

        self.plotbutton.configure(command=plot)
        #self.savebutton.configure(command=save)

    def selecttab(self):
        run = self.gui.ascot.data.active

        # Set values for the coordinate drop-down menus
        quantities = run.getstate_list()
        self.xcrd.setvals(list(quantities.keys()), "r", 0)
        self.ycrd.setvals(list(quantities.keys()), "z", 0)

        # For z and c one can choose "None" which is also the default
        qwithnone = quantities.copy()
        qwithnone["None"] = ("None","")
        self.zcrd.setvals(list(qwithnone.keys()), "None", 0)
        self.ccrd.setvals(list(qwithnone.keys()), "None", 0)

        end, err = run.getstate_markersummary()
        _, end = zip(*end)
        end = ["all"] + list(end)
        self.endc.setvals(end, "all")

        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        """Canvas for plotting the state plots.
        """

        def __init__(self, canvas):
            """Add attribute axes3d which flags whether next axes should be 3d.
            """
            self.axes3d = False
            super().__init__(canvas)

        def make3d(self):
            """Flag next axes to be created 3d.
            """
            self.axes3d = True

        def set_axes(self):
            """Override this method to allow generation of 3d axes.
            """
            if self.axes3d:
                self.axes3d = False
                return self.fig.add_subplot(1,1,1, projection="3d")

            return self.fig.add_subplot(1,1,1)

class Poincare(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)
        self.plot = Poincare.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)

        f1.grid(row=0, column=0, sticky="w")
        f2.grid(row=0, column=1, sticky="ne")

        self.columnconfigure(1, weight=1)

        self.plane = DropdownMenu(f1, label="           Plane: ", width=25,
                                  labelwidth=10, labelanchor="w")
        #self.coord = DropdownMenu(f1, label="Coordinates: ", width=8,
        #                          labelwidth=10, labelanchor="w")
        self.plane.pack(fill="x", anchor="w")
        #self.coord.pack(fill="x", anchor="w")

        self.plotbutton = tk.Button(f2, text="Plot", width=3)
        self.savebutton = tk.Button(f2, text="Store", width=3)
        self.plotbutton.pack(anchor="nw")
        self.savebutton.pack(anchor="nw")

        def plot():
            run = self.gui.ascot.data.active
            self.plot.clear()
            plane = self.plane.getval().split(" ")
            plane = plane[0] + " " + plane[1]

            run.plotorbit_poincare(plane, connlen=True, axes=self.plot.axes)
            self.plot.draw()

        self.plotbutton.configure(command=plot)
        #self.savebutton.configure(command=save)


    def selecttab(self):
        run = self.gui.ascot.data.active
        pol, tor, rad = run.getorbit_poincareplanes()
        planes = []
        for i, p in enumerate(pol):
            planes += ["pol %d at %f deg" % (i+1, p[0])]
        for i, p in enumerate(tor):
            planes += ["tor %d at %f deg" % (i+1, p[0])]
        for i, p in enumerate(rad):
            planes += ["rad %d at %f rho" % (i+1, p[0])]

        self.plane.setvals(planes, planes[0])
        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        pass

class Dists(ContentTab):

    def __init__(self, frame, canvas, gui):
        """Initializes frame by creating all the widgets.
        """
        super().__init__(frame, canvas, gui)
        self.plot = Dists.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)

        f1.grid(row=0, column=0, sticky="w")
        f2.grid(row=0, column=1, sticky="ne")

        self.columnconfigure(1, weight=1)

        self.source = DropdownMenu(f1, label="Source: ", width=24,
                                   labelwidth=8, labelanchor="w",
                                   trace=self.setcoords)
        self.xcrd = DropdownMenu(f1, label="x: ", width=5, labelwidth=2,
                                 trace=self.setxcoord)
        self.ycrd = DropdownMenu(f1, label="y: ", width=5, labelwidth=2)
        self.source.grid(row=0, column=0)
        self.xcrd.grid(row=1, column=0, sticky="e")
        self.ycrd.grid(row=2, column=0, sticky="e")

        self.plotbutton = tk.Button(f2, text="Plot", width=3)
        self.savebutton = tk.Button(f2, text="Store", width=3)
        self.plotbutton.pack(anchor="nw")
        self.savebutton.pack(anchor="nw")

        def plot():
            whichdist = self.source.getval()
            run  = self.gui.ascot.data.active
            dist = run.getdist(whichdist)
            x = self.xcrd.getval()
            y = self.ycrd.getval()

            integrate = {}
            for a in dist.abscissae:
                if x != a and y != a:
                    integrate[a] = np.s_[:]
            dist.integrate(**integrate)

            self.plot.clear()
            run.plotdist(dist, axes=self.plot.axes)
            self.plot.draw()

        self.plotbutton.configure(command=plot)

    def setcoords(self, *args):
        dist = self.source.var.get()
        coords = []
        if "com" in dist:
            coords += ["mu", "ekin", "ptor"]
        else:
            if "rho" in dist:
                coords += ["rho", "theta", "phi"]
            else:
                coords += ["r", "z", "phi"]

            if "5d" in dist:
                coords += ["ppar", "pperp"]
            else:
                coords += ["pr", "pz", "pphi"]

            coords += ["time", "charge"]

        self.xcrd.setvals(coords, coords[0])
        self.ycrd.setvals(coords[1:]+["None"], "None")

    def setxcoord(self, *args):
        dist = self.source.var.get()
        x = self.xcrd.var.get()
        y = self.ycrd.var.get()

        coords = []
        if "com" in dist:
            coords += ["mu", "ekin", "ptor"]
        else:
            if "rho" in dist:
                coords += ["rho", "theta", "phi"]
            else:
                coords += ["r", "z", "phi"]

            if "5d" in dist:
                coords += ["ppar", "pperp"]
            else:
                coords += ["pr", "pz", "pphi"]

            coords += ["time", "charge"]

        coords.remove(x)
        if x == y:
            self.ycrd.setvals(coords, "None")
        else:
            self.ycrd.setvals(coords, y)

    def selecttab(self):
        run = self.gui.ascot.data.active
        dists = run.getdist(None)
        self.source.setvals(dists, dists[0])

        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        pass

class Moments(ContentTab):

    def __init__(self, frame, canvas, gui):
        """Initializes frame by creating all the widgets.
        """
        super().__init__(frame, canvas, gui)
        self.plot = Moments.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)
        f1.grid(row=0, column=0, sticky="w")
        f2.grid(row=0, column=1, sticky="ne")

        self.columnconfigure(1, weight=1)

        self.source = DropdownMenu(f1, label="  Source: ", width=24,
                                   labelwidth=8, labelanchor="w")
        self.qnt = DropdownMenu(f1, label="Quantity: ", width=24,
                                labelwidth=8, labelanchor="w")
        self.source.grid(row=0, column=0)
        self.qnt.grid(row=1, column=0)

        self.plotbutton = tk.Button(f2, text="Plot", width=3)
        self.savebutton = tk.Button(f2, text="Store", width=3)
        self.plotbutton.pack(anchor="nw")
        self.savebutton.pack(anchor="nw")

        self.qnt.setvals(["density"], "density")

        def plot():
            dist = self.source.var.get()
            moment  = self.qnt.var.get()
            run  = self.gui.ascot.data.active
            dist = run.getdist(dist)
            mom  = run.getdist_moments(dist, moment)

            self.plot.clear()
            run.plotdist_moments(mom, moment, axes=self.plot.axes)
            self.plot.draw()

        self.plotbutton.configure(command=plot)

    def selecttab(self):
        run = self.gui.ascot.data.active
        dists = run.getdist(None)
        vals = []
        if "5d" in dists: vals += ["5d"]
        if "rho5d" in dists: vals += ["rho5d"]
        self.source.setvals(vals, vals[0])

        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        pass

class LossSummary(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)
        self.plot = LossSummary.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        self.losssummary = tk.Text(self, width=50)
        self.losssummary.pack(fill="both", expand=True)
        self.losssummary.config(state="disabled")

    def selecttab(self):
        run = self.gui.ascot.data.active
        self.losssummary.config(state="normal")
        self.losssummary.delete("1.0", "end")

        text = ""
        msg = run.getstate_losssummary()
        for m in msg:
            text += m + "\n"

        self.losssummary.insert("end", text)
        self.losssummary.config(state="disabled")

        self.plot.clear()
        #run.endstate.plot_lossmap(5.0, axes=self.canvas.axes[4])
        self.plot.draw()
        self.canvas.slideshow(self)

    class Canvas(PlotFrame):

        def set_axes(self):
            """Override this method to allow generation of multiple axes.
            """
            gs = GridSpec(2,3)
            ax = []
            ax += [self.fig.add_subplot(gs[:,0])]
            ax += [self.fig.add_subplot(gs[0,1])]
            ax += [self.fig.add_subplot(gs[1,1])]
            ax += [self.fig.add_subplot(gs[0,2])]
            ax += [self.fig.add_subplot(gs[1,2])]
            return ax

class WallLoad(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)
        self.plot = WallLoad.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        self.loadsummary = tk.Text(self, width=50)
        self.loadsummary.pack(fill="both", expand=True)
        self.loadsummary.config(state="disabled")

    def selecttab(self):
        run = self.gui.ascot.data.active
        self.loadsummary.config(state="normal")
        self.loadsummary.delete("1.0", "end")

        wetted, epeak = run.getwall_figuresofmerit()
        text = "Wetted area: " \
            + np.format_float_scientific(wetted, precision=2) \
            + " " + str(wetted.units) + "\n"\
            + "Peak load: " \
            + np.format_float_scientific(epeak, precision=2) \
            + " " + str(epeak.units) + "\n"\

        self.loadsummary.insert("end", text)
        self.loadsummary.config(state="disabled")

        self.plot.clear()
        run.plotwall_loadvsarea(axes=self.plot.axes[0])
        run.plotwall_torpol(axes=self.plot.axes[1])
        self.plot.draw()
        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        """Canvas for plotting.
        """

        def set_axes(self):
            """Override this method to allow generation of multiple axes.
            """
            gs = GridSpec(1,2)
            ax = []
            ax += [self.fig.add_subplot(gs[:,0])]
            ax += [self.fig.add_subplot(gs[0,1])]
            return ax


class Wall3D(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)
        self.plot = Wall3D.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        # Use three frames for layout
        f1 = ttk.Frame(self) # Camera controls
        f2 = ttk.Frame(self) # Plot and save buttons
        f3 = ttk.Frame(f2) # Other settings (located inside f2)

        f1.grid(row=0, column=0, sticky="w")
        f2.grid(row=0, column=1, sticky="nwes")

        self.columnconfigure(1, weight=1)

        # Camera selection
        self.csel = DropdownMenu(f1, label="", width=12, labelwidth=0)
        self.csel.grid(row=0, column=0, columnspan=2, sticky="w")
        self.csel.setvals(["Camera 1"],#, "Camera 2", "Camera 3",
                          #"Camera 4", "Camera 5", "Camera 6",
                          #"Camera 7", "Camera 8", "Camera 9"],
                          "Camera 1")

        # Interactive plot
        self.intrbutton = tk.Button(f1, text="Interactive", width=6)
        self.intrbutton.grid(row=0, column=2, sticky="e")

        # Camera position control
        l = tk.Label(f1, text="Camera location:")
        l.grid(row=1, column=0, columnspan=2, sticky="w")
        self.cposr = NumEntry(f1, labeltext="R", entrywidth=7,
                              labelwidth=2, anchor="e", defval=0.0)
        self.cposp = NumEntry(f1, labeltext=" phi", entrywidth=7,
                              labelwidth=4, anchor="e", defval=0.0)
        self.cposz = NumEntry(f1, labeltext=" z", entrywidth=7,
                              labelwidth=2, anchor="e", defval=0.0)

        self.cposr.grid(row=2, column=0)
        self.cposp.grid(row=2, column=1)
        self.cposz.grid(row=2, column=2)

        # Camera focus control
        l = tk.Label(f1, text="Focus point:")
        l.grid(row=3, column=0, columnspan=2, sticky="w")
        self.cfocr = NumEntry(f1, labeltext="R", entrywidth=7,
                              labelwidth=2, anchor="e", defval=0.0)
        self.cfocp = NumEntry(f1, labeltext="phi", entrywidth=7,
                              labelwidth=4, anchor="e", defval=0.0)
        self.cfocz = NumEntry(f1, labeltext=" z", entrywidth=7,
                              labelwidth=2, anchor="e", defval=0.0)

        self.cfocr.grid(row=4, column=0)
        self.cfocp.grid(row=4, column=1)
        self.cfocz.grid(row=4, column=2)

        # Camera angle control
        l = tk.Label(f1, text="Azimuth:")
        l.grid(row=5, column=0, sticky="w")
        l = tk.Label(f1, text="Elevation:")
        l.grid(row=5, column=1, sticky="w")
        l = tk.Label(f1, text="Rotation:")
        l.grid(row=5, column=2, sticky="w")
        self.canga = NumEntry(f1, labeltext=" ", entrywidth=7,
                              labelwidth=2, anchor="e", defval=0.0)
        self.cange = NumEntry(f1, labeltext="   ", entrywidth=7,
                              labelwidth=4, anchor="e", defval=0.0)
        self.cangr = NumEntry(f1, labeltext="  ", entrywidth=7,
                              labelwidth=2, anchor="e", defval=0.0)

        self.canga.grid(row=6, column=0)
        self.cange.grid(row=6, column=1)
        self.cangr.grid(row=6, column=2)

        # Plot, save, and interactive buttons
        self.plotbutton = tk.Button(f2, text="Plot", width=3)
        self.savebutton = tk.Button(f2, text="Store", width=3)
        self.plotbutton.pack(anchor="ne")
        self.savebutton.pack(anchor="ne")

        f3.pack(fill="both")
        l = tk.Label(f3, text="Markers")
        l.grid(row=0, column=1, sticky="w")
        l = tk.Label(f3, text="Wall loads")
        l.grid(row=1, column=1, sticky="w")
        l = tk.Label(f3, text="log10 scale")
        l.grid(row=2, column=1, sticky="w")
        self.showmarkers = Tickbox(f3, width=1)
        self.showmarkers.grid(row=0, column=0)
        self.showloads = Tickbox(f3, width=1)
        self.showloads.grid(row=1, column=0)
        self.logscale = Tickbox(f3, width=1)
        self.logscale.grid(row=2, column=0)

        self.wallmesh = None

    def selecttab(self):
        """Setup data.
        """
        run = self.gui.ascot.data.active
        self.wallmesh = run.getwall_3dmesh()
        self.setcamera(*defaultcamera(self.wallmesh))

        self.plotbutton.configure(command=lambda:self.still(run))
        self.intrbutton.configure(command=lambda:self.interactive(run))

        self.canvas.slideshow(self)

    def still(self, run):
        """Plot still picture on the canvas.
        """
        self.plot.clear()
        cpos, cfoc, cang = self.getcamera()
        points = None
        if self.showmarkers.getval() == 1:
            points = run.getstate_pointcloud(endcond="wall")
        data = "eload" if self.showloads.getval() == 1 else None
        log = self.logscale.getval() == 1
        run.plotwall_3dstill(self.wallmesh,
                             points=points, data=data, log=log,
                             cpos=cpos, cfoc=cfoc, cang=cang,
                             axes=self.plot.axes)
        self.plot.draw()

    def interactive(self, run):
        """Open new window for interactive 3D view.
        """
        def recordcamera(plotter):
            cpos = plotter.camera_position[0]
            cfoc = plotter.camera_position[1]
            cang = plotter.camera_position[2]
            cang = [plotter.camera.azimuth,
                    plotter.camera.elevation,
                    plotter.camera.roll]

            self.setcamera(cpos, cfoc, cang)

        cpos, cfoc, cang = self.getcamera()
        points = None
        if self.showmarkers.getval() == 1:
            points = run.getstate_pointcloud(endcond="wall")
        data = "eload" if self.showloads.getval() == 1 else None
        log = self.logscale.getval() == 1
        run.plotwall_3dinteractive(self.wallmesh,
                                   ("k", recordcamera),
                                   points=points, data=data, log=log,
                                   cpos=cpos, cfoc=cfoc, cang=cang)

    def setcamera(self, cpos, cfoc, cang):
        """Record and display camera coordinates (given in xyz) in rpz.
        """

        def roundval(val):
            """
            Round to three decimals.
            """
            return int(val * 100) / 100

        cr = np.sqrt( cpos[0]**2 + cpos[1]**2 )
        cp = np.arctan2( cpos[1], cpos[0] )
        cz = cpos[2]

        self.cposr.setval( roundval(cr) )
        self.cposp.setval( roundval(cp * 180 / np.pi) )
        self.cposz.setval( roundval(cz) )

        cr = np.sqrt( cfoc[0]**2 + cfoc[1]**2 )
        cp = np.arctan2( cfoc[1], cfoc[0] )
        cz = cfoc[2]

        self.cfocr.setval( roundval(cr) )
        self.cfocp.setval( roundval(cp * 180 / np.pi) )
        self.cfocz.setval( roundval(cz) )

        self.canga.setval( roundval(cang[0]) )
        self.cange.setval( roundval(cang[1]) )
        self.cangr.setval( roundval(cang[2]) )

    def getcamera(self):
        """Get camera coordinates (stored as rpz) as xyz.
        """
        cr = self.cposr.getval()
        cp = self.cposp.getval() * np.pi / 180
        cz = self.cposz.getval()

        cpos = [cr * np.cos(cp), cr * np.sin(cp), cz]

        cr = self.cfocr.getval()
        cp = self.cfocp.getval() * np.pi / 180
        cz = self.cfocz.getval()

        cfoc = [cr * np.cos(cp), cr * np.sin(cp), cz]

        cang = [self.canga.getval(),
                self.cange.getval(),
                self.cangr.getval()]

        return (cpos, cfoc, cang)

    class Canvas(PlotFrame):
        pass
