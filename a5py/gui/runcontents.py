import tkinter as tk
from tkinter import ttk

from .components import ContentTab, PlotFrame, NumEntry, Tickbox, DropdownMenu

from a5py.ascot5io.options import Opt
from a5py.ascot5io.marker import Marker

class Trace(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)
        self.plot = Trace.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)
        f2.pack(side="left", fill="both", expand=True)
        f1.pack(side="right", fill="y")


        # Right-hand side
        f1.columnconfigure(1, weight=1)

        self.plotbutton = tk.Button(f1, text="Run and plot", width=12)
        self.plotbutton.grid(column=1, row=0, rowspan=2, sticky="e")
        self.gamemode = DropdownMenu(f1, label="", width=15,
                                     labelwidth=0)
        self.gamemode.grid(column=1, row=2, rowspan=2, sticky="e")
        self.gamemode.setvals(["Guding centre", "Full gyro-orbit"],
                              "Guiding centre")
        self.axeseq = Tickbox(f1, label=" Axes equal", width=12)
        self.axeseq.grid(column=1, row=4, sticky="e")

        self.milent = NumEntry(f1, labeltext="Mileage [s]", entrywidth=10,
                               labelwidth=12, anchor="w", defval=0.00001)
        self.pntent = NumEntry(f1, labeltext="N points", entrywidth=10,
                               labelwidth=12, anchor="w", defval=1000,
                               isint=True)

        self.milent.grid(column=1, row=6)
        self.pntent.grid(column=1, row=7)

        # Left-hand side
        self.species = DropdownMenu(f2, label=" species: ", width=8,
                                     labelwidth=10, num_dropdown_items=15)
        self.species.grid(column=0, row=0, rowspan=2, sticky="e")
        self.species.setvals(["e", "H", "D", "T", "He3", "He4", "Be9", "C12",
                               "Ne20", "Ar40", "Ni59", "Xe132", "W184"], "e")

        self.r_entry = NumEntry(f2, labeltext=" R [m]", entrywidth=10,
                                labelwidth=14, anchor="e", defval=6.7)
        self.p_entry = NumEntry(f2, labeltext=" phi [deg]", entrywidth=10,
                                labelwidth=14, anchor="e", defval=0.0)
        self.z_entry = NumEntry(f2, labeltext=" z [m]", entrywidth=10,
                                labelwidth=14, anchor="e", defval=0.0)
        self.e_entry = NumEntry(f2, labeltext=" Energy [eV]", entrywidth=10,
                                labelwidth=14, anchor="e", defval=10.0e6)
        self.k_entry = NumEntry(f2, labeltext=" Pitch [vpa/v]", entrywidth=10,
                                labelwidth=14, anchor="e", defval=0.99)
        self.r_entry.grid(column=0, row=3)
        self.p_entry.grid(column=0, row=4)
        self.z_entry.grid(column=0, row=5)
        self.e_entry.grid(column=0, row=6)
        self.k_entry.grid(column=0, row=7)

        def plot():
            r      = self.r_entry.getval()
            phi    = self.p_entry.getval()
            z      = self.z_entry.getval()
            energy = self.e_entry.getval()
            pitch  = self.k_entry.getval()
            gamemode = self.gamemode.getval()
            if gamemode == "Guiding centre":
                print("gamemode GC")
                simmode = 2
            elif gamemode == "Full gyro-orbit":
                print("gamemode FO")
                simmode = 1

            mrk = Marker.generate("gc", n=1, species=self.species.getval())
            mrk["r"][:]      = r
            mrk["phi"][:]    = phi
            mrk["z"][:]      = z
            mrk["energy"][:] = energy
            mrk["pitch"][:]  = pitch

            mil = self.milent.getval()
            npoint = self.pntent.getval()
            interval = mil/npoint

            ids = mrk["ids"]

            opt = Opt.get_default()
            opt.update({"SIM_MODE"               : simmode,
                        "ENABLE_ADAPTIVE"        : 1,
                        "ENABLE_ORBIT_FOLLOWING" : 1,
                        "ENDCOND_SIMTIMELIM"     : 1,
                        "ENDCOND_LIM_SIMTIME"    : 100,
                        "ENDCOND_MAX_MILEAGE"    : mil,
                        "ENABLE_ORBITWRITE"      : 1,
                        "ORBITWRITE_MODE"        : 1,
                        "ORBITWRITE_NPOINT"      : npoint,
                        "ORBITWRITE_INTERVAL"    : interval
            })

            self.gui.ascot.simulation_initmarkers(**mrk)
            self.gui.ascot.simulation_initoptions(**opt)
            vrun = self.gui.ascot.simulation_run()

            self.plot.clear()
            vrun.plotorbit_trajectory("r",
                                      "z",
                                      axes=self.plot.axes,
                                      axesequal=self.axeseq.getval(),
                                      )
            self.plot.draw()

            self.gui.ascot.simulation_free(diagnostics=True)

        self.plotbutton.configure(command=plot)

    def selecttab(self):
        self.canvas.slideshow(self)

    class Canvas(PlotFrame):
        pass
