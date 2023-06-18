import tkinter as tk
from tkinter import ttk

from collections import OrderedDict
from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class OrbitFrame(ttk.Frame):

    def init(self, gui, canvas):
        self.gui = gui

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
            ("|v|",       "m/s"),
            ("vpa",       "m/s"),
            ("vpe",       "m/s"),
            ("vR",        "m/s"),
            ("vphi",      "m/s"),
            ("vz",        "m/s"),
            ("|B|",       "T"),
            ("BR",        "T"),
            ("Bphi",      "T"),
            ("Bz",        "T"),
            ("charge",    "T")
        ])

        class TrajectoryFrame(ttk.Frame):

            def init(self, canvas):
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

                self.xcrd.setvals(list(quantities.keys()), "R", 0)
                self.ycrd.setvals(list(quantities.keys()), "z", 0)
                qwithnone = quantities.copy()
                qwithnone["None"] = ("None","")
                self.zcrd.setvals(list(qwithnone.keys()), "None", 0)
                self.ccrd.setvals(list(qwithnone.keys()), "None", 0)

                self.canvas = canvas
                return self


            def plot(self, run):
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

                if zcoord == "None":
                    self.canvas.fig_rzview.clear()

                    if ccoord == "None":
                        run.orbit.plot(xcoord, ycoord, equal=equal, log=log,
                                          endcond=endcond, axes=self.canvas.fig_rzview.axes)
                    else:
                        run.orbit.scatter(xcoord, ycoord, c=ccoord, equal=equal,
                                             log=log, endcond=endcond, axes=self.canvas.fig_rzview.axes)
                else:
                    self.canvas.fig_rzview.make3d()
                    self.canvas.fig_rzview.clear()
                    if ccoord == "None":
                        run.orbit.plot(xcoord, ycoord, zcoord, equal=equal,
                                          log=log, endcond=endcond, axes=self.canvas.fig_rzview.axes)
                    else:
                        run.orbit.scatter(xcoord, ycoord, zcoord, ccoord,
                                             equal=equal, log=log,
                                             endcond=endcond, axes=self.canvas.fig_rzview.axes)

                self.canvas.fig_rzview.draw()


        class PoincareFrame(ttk.Frame):

            def init(self, canvas):

                f1 = ttk.Frame(self)
                f2 = ttk.Frame(self)

                f1.grid(row=0, column=0, sticky="w")
                f2.grid(row=0, column=1, sticky="ne")

                self.columnconfigure(1, weight=1)

                self.plane = DropdownMenu(f1, label="           Plane: ", width=8,
                                          labelwidth=10, labelanchor="w")
                #self.coord = DropdownMenu(f1, label="Coordinates: ", width=8,
                #                          labelwidth=10, labelanchor="w")
                self.plane.pack(fill="x", anchor="w")
                #self.coord.pack(fill="x", anchor="w")

                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                self.canvas = canvas
                return self

            def plot(self, run):
                self.canvas.fig_rzview.clear()
                plane = self.plane.getval()
                plane = run.getorbit_poincareplanes().index(plane)

                run.plotorbit_poincare(plane, conlen=True, axes=self.canvas.fig_rzview.axes)
                self.canvas.fig_rzview.draw()


        master = ttk.Notebook(self)

        self.frametrajectory = TrajectoryFrame(master).init(canvas)
        self.framepoincare   = PoincareFrame(master).init(canvas)
        master.add(self.frametrajectory, text="Trajectory")
        master.add(self.framepoincare,   text="Poincare")
        master.pack(fill="both", expand=True)

        self.canvas = canvas

        return self


    def display(self):
        run = self.gui.ascot.data.active
        try:
            endconds, counts = run.endstate.listendconds()
        except AttributeError:
            endconds, counts = run.inistate.listendconds()

        endconds = ["all"] + endconds

        self.frametrajectory.endc.setvals(endconds, "all")
        planes = run.getorbit_poincareplanes()
        self.framepoincare.plane.setvals(planes, planes[0])

        def plottrajectory():
            self.frametrajectory.plot(run)

        def plotpoincare():
            self.framepoincare.plot(run)

        self.frametrajectory.plotbutton.configure(command=plottrajectory)
        self.framepoincare.plotbutton.configure(command=plotpoincare)


class OrbitCanvas(ttk.Frame):

    def init(self):
        class Plot3DFrame(PlotFrame):

            def __init__(self, master, tight_layout=True):
                self.axes3d = False
                super().__init__(master, tight_layout=tight_layout)

            def make3d(self):
                self.axes3d = True

            def set_axes(self):
                if self.axes3d:
                    self.axes3d = False
                    return self.fig.add_subplot(1,1,1, projection="3d")

                return self.fig.add_subplot(1,1,1)

        fig_rzview = Plot3DFrame(self)
        fig_rzview.place(relheight=0.8, anchor="nw")

        self.fig_rzview = fig_rzview
        return self
