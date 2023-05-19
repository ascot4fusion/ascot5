import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry

from a5py.ascot5io import options
from a5py.ascot5io import mrk as marker

class OrbitFrame(ttk.Frame):

    def init(self, canvas):
        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)
        f1.pack(fill="x")
        f2.pack(fill="both")

        f1.columnconfigure(1, weight=1)

        self.milent = NumEntry(f1, labeltext="Mileage [s]", entrywidth=10,
                               labelwidth=12, anchor="w", defval=0.0001)
        self.pntent = NumEntry(f1, labeltext="N points", entrywidth=10,
                               labelwidth=12, anchor="w", defval=1000,
                               isint=True)

        self.milent.grid(column=0, row=0)
        self.pntent.grid(column=0, row=1)

        self.plotbutton = tk.Button(f1, text="Run and plot", width=8)
        self.plotbutton.grid(column=1, row=0, rowspan=2)

        self.r_entry = NumEntry(f2, labeltext="R [m]", entrywidth=5,
                                labelwidth=12, anchor="w", defval=6.7)
        self.p_entry = NumEntry(f2, labeltext="phi [deg]", entrywidth=5,
                                labelwidth=12, anchor="w", defval=0.0)
        self.z_entry = NumEntry(f2, labeltext="z [m]", entrywidth=5,
                                labelwidth=12, anchor="w", defval=0.0)
        self.e_entry = NumEntry(f2, labeltext="Energy [eV]", entrywidth=10,
                                labelwidth=14, anchor="w", defval=10.0e6)
        self.k_entry = NumEntry(f2, labeltext="Pitch [vpa/v]", entrywidth=10,
                                labelwidth=14, anchor="w", defval=0.99)
        self.r_entry.grid(column=0, row=0)
        self.p_entry.grid(column=0, row=1)
        self.z_entry.grid(column=0, row=2)
        self.e_entry.grid(column=1, row=0)
        self.k_entry.grid(column=1, row=1)

        self.canvas     = canvas


    def display(self, gui):

        gui.ascot.simulation_initinput()

        def plot():
            r      = self.r_entry.getval()
            phi    = self.p_entry.getval()
            z      = self.z_entry.getval()
            energy = self.e_entry.getval()
            pitch  = self.k_entry.getval()

            mrk = marker.generatemrk(nmrk=1, mrktype="gc", species="electron")
            mrk["r"][:]      = r
            mrk["phi"][:]    = phi
            mrk["z"][:]      = z
            mrk["energy"][:] = energy
            mrk["pitch"][:]  = pitch

            mil = 1e-5
            npoint = 1000
            interval = mil/npoint

            ids = mrk["ids"]

            opt = options.generateopt(clean=True)
            opt.update({"SIM_MODE"               : 2,
                        "ENABLE_ADAPTIVE"        : 1,
                        "ENABLE_ORBIT_FOLLOWING" : 1,
                        "ENDCOND_SIMTIMELIM"     : 1,
                        "ENDCOND_MAX_SIMTIME"    : 100,
                        "ENDCOND_MAX_MILEAGE"    : mil,
                        "ENABLE_ORBITWRITE"      : 1,
                        "ORBITWRITE_MODE"        : 1,
                        "ORBITWRITE_NPOINT"      : npoint,
                        "ORBITWRITE_INTERVAL"    : interval
            })

            gui.ascot.simulation_initmarkers("gc", mrk)
            gui.ascot.simulation_initoptions(opt)
            gui.ascot.simulation_run()
            rorb = gui.ascot.simulation_getorbit("r", ids)
            zorb = gui.ascot.simulation_getorbit("z", ids)
            gui.ascot.simulation_freeoutput()

            self.canvas.fig_rzview.clear()
            self.canvas.fig_rzview.axes.plot(rorb, zorb)
            self.canvas.fig_rzview.draw()

        self.plotbutton.configure(command=plot)


class OrbitCanvas(ttk.Frame):

    def init(self):
        fig_rzview = PlotFrame(self)
        fig_rzview.place(relheight=0.8, anchor="nw")

        self.fig_rzview = fig_rzview
        return self
