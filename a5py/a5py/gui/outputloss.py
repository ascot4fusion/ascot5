import tkinter as tk
from tkinter import ttk
import numpy as np
import pyvista as pv

from matplotlib.gridspec import GridSpec
from a5py.plotting import defaultcamera
from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class LossFrame(ttk.Notebook):

    def init(self, gui, canvas):
        self.gui = gui

        class GeneralFrame(ttk.Frame):

            def init(self, canvas):
                self.losssummary = tk.Text(self, width=50)
                self.losssummary.pack(fill="both", expand=True)
                self.losssummary.config(state="disabled")

                self.run = None
                self.canvas = canvas
                return self


            def load(self, run):
                self.run = run
                self.losssummary.config(state="normal")
                self.losssummary.delete("1.0", "end")

                text = ""
                msg = run.getstate_losssummary()
                for m in msg:
                    text += m + "\n"

                self.losssummary.insert("end", text)
                self.losssummary.config(state="disabled")


            def plot(self):
                self.canvas.prepare("general")
                self.canvas.clear()
                if self.run is not None:
                    self.run.endstate.plot_lossmap(5.0, axes=self.canvas.axes[4])
                self.canvas.draw()


        class WallloadFrame(ttk.Frame):

            def init(self, canvas):
                self.run = None
                self.loadsummary = tk.Text(self, width=50)
                self.loadsummary.pack(fill="both", expand=True)
                self.loadsummary.config(state="disabled")

                self.canvas = canvas
                return self


            def load(self, run):
                self.run = run
                self.loadsummary.config(state="normal")
                self.loadsummary.delete("1.0", "end")

                text = ""
                msg = run.getwall_figuresofmerit()
                for m in msg:
                    text += m + "\n"

                self.loadsummary.insert("end", text)
                self.loadsummary.config(state="disabled")


            def plot(self):
                self.canvas.prepare("wallload")
                self.canvas.clear()
                self.run.plotwall_loadvsarea(axes=self.canvas.axes[0])
                self.canvas.draw()


        class View3DFrame(ttk.Frame):

            def init(self, canvas):
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
                self.points   = None
                self.canvas   = canvas

                return self


            def load(self, run):
                """
                Setup data.
                """

                self.wallmesh = run.getwall_3dmesh()
                self.points = run.getstate_pointcloud(endcond="wall")

                self.setcamera(*defaultcamera(self.wallmesh))

                self.plotbutton.configure(command=lambda:self.still(run))
                self.intrbutton.configure(command=lambda:self.interactive(run))


            def still(self, run):
                """
                Plot still picture on the canvas.
                """
                self.canvas.prepare("")
                self.canvas.clear()
                cpos, cfoc, cang = self.getcamera()
                points = self.points if self.showmarkers.getval() == 1 else None
                data = "eload" if self.showloads.getval() == 1 else None
                log = self.logscale.getval() == 1
                run.plotwall_3dstill(self.wallmesh,
                                     points=points, data=data, log=log,
                                     cpos=cpos, cfoc=cfoc, cang=cang,
                                     axes=self.canvas.axes[0])
                self.canvas.draw()


            def interactive(self, run):
                """
                Open new window for interactive 3D view.
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
                points = self.points if self.showmarkers.getval() == 1 else None
                data = "eload" if self.showloads.getval() == 1 else None
                log = self.logscale.getval() == 1
                run.plotwall_3dinteractive(self.wallmesh,
                                           ("k", recordcamera),
                                           points=points, data=data, log=log,
                                           cpos=cpos, cfoc=cfoc, cang=cang)


            def setcamera(self, cpos, cfoc, cang):
                """
                Record and display camera coordinates (given in xyz) in rpz.
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
                """
                Get camera coordinates (stored as rpz) as xyz.
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


        #master = ttk.Notebook(self)
        self.framegeneral  = GeneralFrame(self).init(canvas)
        self.framewallload = WallloadFrame(self).init(canvas)
        self.frameview3d   = View3DFrame(self).init(canvas)
        self.add(self.framegeneral,  text="General")
        self.add(self.framewallload, text="Wall loads")
        self.add(self.frameview3d,   text="View 3D")
        #self.pack(fill="both", expand=True)

        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            if tab == "General":
                self.framegeneral.plot()
            elif tab == "Wall loads":
                self.framewallload.plot()
            elif tab == "View 3D":
                pass

        self.bind('<<NotebookTabChanged>>', on_tab_change)

        self.canvas = canvas
        return self


    def display(self):
        run = self.gui.ascot.hdf5.active
        self.framegeneral.load(run)
        self.framewallload.load(run)
        self.frameview3d.load(run)


class LossCanvas(PlotFrame):

    def __init__(self, master):
        """
        Add attribute.
        """
        self.frame = ""
        super().__init__(master)

    def init(self):
        self.general_phitheta = 0
        self.general_energy   = 1
        self.general_mileage  = 2
        self.general_lossmap1 = 3
        self.general_lossmap2 = 4
        return self


    def prepare(self, frame):
        self.frame = frame


    def set_axes(self):
        """
        Override the PlotFrame method in order to create multiple axes.
        """
        ax = []
        if self.frame == "general":
            gs = GridSpec(2,3)
            ax += [self.fig.add_subplot(gs[:,0])]
            ax += [self.fig.add_subplot(gs[0,1])]
            ax += [self.fig.add_subplot(gs[1,1])]
            ax += [self.fig.add_subplot(gs[0,2])]
            ax += [self.fig.add_subplot(gs[1,2])]
        elif self.frame == "wallload":
            ax += [self.fig.add_subplot(1,1,1)]
        elif self.frame == "":
            ax += [self.fig.add_subplot(1,1,1)]

        return ax
