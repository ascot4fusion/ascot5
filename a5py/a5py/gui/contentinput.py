import tkinter as tk
from tkinter import ttk
import numpy as np

from .components import PlotFrame, NumEntry
#from.guiparams import retrieve, store
from a5py.ascotpy.libbfield  import LibBfield
from a5py.ascotpy.libefield  import LibEfield
from a5py.ascotpy.libplasma  import LibPlasma
from a5py.ascotpy.libneutral import LibNeutral
from a5py.ascotpy.libboozer  import LibBoozer
from a5py.ascotpy.libmhd     import LibMhd

class ContentInput:

    def __init__(self, gui, settings, canvas):
        self.gui = gui

        ## Add widgets

        class EntryPanel(tk.Frame):

            def init(self):
                xmin_entry = NumEntry(panel, labeltext="R = ", defval=0.1)
                xmax_entry = NumEntry(panel, labeltext=" - ", defval=10.0)
                xnum_entry = NumEntry(panel, labeltext=" : ", defval=50,
                                      isint=True)

                ymin_entry = NumEntry(panel, labeltext="z = ", defval=8.0)
                ymax_entry = NumEntry(panel, labeltext=" - ", defval=-8.0)
                ynum_entry = NumEntry(panel, labeltext=" : ", defval=50,
                                      isint=True)

                xmin_entry.grid(row=0, column=0, sticky="W")
                xmax_entry.grid(row=0, column=1, sticky="W")
                xnum_entry.grid(row=0, column=2, sticky="W")

                ymin_entry.grid(row=1, column=0, sticky="W")
                ymax_entry.grid(row=1, column=1, sticky="W")
                ynum_entry.grid(row=1, column=2, sticky="W")

                phi_entry = NumEntry(panel, labeltext="phi [deg]:", defval=0)
                phi_entry.grid(row=3, column=0)

                time_entry = NumEntry(panel, labeltext="time [s]:", defval=0)
                time_entry.grid(row=3, column=1)

                cmin_entry = NumEntry(panel, labeltext="Color range", defval="")
                cmax_entry = NumEntry(panel, labeltext=" - ", defval="")
                cmin_entry.grid(row=2, column=0)
                cmax_entry.grid(row=2, column=1)

                qchoice = tk.StringVar(panel)
                qinput = ttk.Combobox(panel, width=10, textvariable=qchoice)
                qinput.grid(row=4, column=1, sticky="WE")

                self.xmin_entry = xmin_entry
                self.xmax_entry = xmax_entry
                self.xnum_entry = xnum_entry
                self.ymin_entry = ymin_entry
                self.ymax_entry = ymax_entry
                self.ynum_entry = ynum_entry
                self.phi_entry  = phi_entry
                self.time_entry = time_entry
                self.cmin_entry = cmin_entry
                self.cmax_entry = cmax_entry
                self.qchoice    = qchoice
                self.qinput     = qinput


            def plotparams(self):
                out = {}
                out["r"] = np.linspace( self.xmin_entry.getval() ,
                                        self.xmax_entry.getval() ,
                                        self.xnum_entry.getval() )

                out["z"] = np.linspace( self.ymin_entry.getval() ,
                                        self.ymax_entry.getval() ,
                                        self.ynum_entry.getval() )

                out["phi"]  = self.phi_entry.getval() * np.pi / 180
                out["time"] = self.time_entry.getval()

                clim = [None, None]
                if not self.cmin_entry.isempty():
                    clim[0] = self.cmin_entry.getval()
                if not self.cmax_entry.isempty():
                    clim[1] = self.cmax_entry.getval()
                out["clim"] = clim
                out["qnt"]  = self.qchoice.get()
                return out


        settingsframe = tk.Frame(settings)
        panel = EntryPanel(settingsframe)
        panel.init()

        savebutton = tk.Button(panel, text="Store")
        savebutton.grid(row=2, column=2, sticky="WE")
        plotbutton = tk.Button(panel, text="Plot")
        plotbutton.grid(row=3, column=2, sticky="WE")

        # Store parameters
        gui.params.add(
            input_rzplot_minr=panel.xmin_entry.choice,
            input_rzplot_maxr=panel.xmax_entry.choice,
            input_rzplot_numr=panel.xnum_entry.choice,
            input_rzplot_minz=panel.ymin_entry.choice,
            input_rzplot_maxz=panel.ymax_entry.choice,
            input_rzplot_numz=panel.ynum_entry.choice,
            input_rzplot_qnt=panel.qchoice
        )

        plotcanvas = tk.Frame(canvas)
        fig_rzview = PlotFrame(plotcanvas)
        fig_rzview.place(relheight=0.8, anchor="nw")

        self.settingsframe = settingsframe
        self.panel         = panel
        self.savebutton    = savebutton
        self.plotbutton    = plotbutton
        self.plotcanvas    = plotcanvas
        self.fig_rzview    = fig_rzview


    def display(self):

        self.settingsframe.pack_forget()
        self.panel.pack_forget()
        self.plotcanvas.pack_forget()

        # Find what quantities we can plot and init
        quantities = []

        bfield = False
        if hasattr(self.gui.ascot.hdf5, "bfield"):
            quantities += LibBfield.quantities
            bfield = True
        efield = False
        if hasattr(self.gui.ascot.hdf5, "efield"):
            quantities += LibEfield.quantities
            efield = True
        plasma = False
        if hasattr(self.gui.ascot.hdf5, "plasma"):
            quantities += LibPlasma.quantities
            plasma = True
        neutral = False
        if hasattr(self.gui.ascot.hdf5, "neutral"):
            quantities += LibNeutral.quantities
            neutral = True
        boozer = False
        if hasattr(self.gui.ascot.hdf5, "boozer"):
            quantities += LibBoozer.quantities
            boozer = True
        mhd = False
        if hasattr(self.gui.ascot.hdf5, "mhd"):
            quantities += LibMhd.quantities
            mhd = True

        self.gui.ascot.init(bfield=bfield, efield=efield, neutral=neutral,
                            plasma=plasma, boozer=boozer, mhd=mhd,
                            ignorewarnings=True)

        self.panel.qinput["values"] = quantities
        if self.panel.qchoice.get() is None or \
           ( self.panel.qchoice.get() not in quantities ):
            self.panel.qchoice.set(quantities[0])

        def storesettings():
            self.gui.params.store(self.gui.ascot.h5fn,
                ["input_rzplot_minr",
                 "input_rzplot_maxr",
                 "input_rzplot_numr",
                 "input_rzplot_minz",
                 "input_rzplot_maxz",
                 "input_rzplot_numz",
                 "input_rzplot_qnt"]
            )

        def plot():
            self.fig_rzview.clear()
            out = self.panel.plotparams()
            self.gui.ascot.plotRz(out["r"], out["phi"], out["z"], out["time"],
                                  out["qnt"], clim=out["clim"],
                                  axes=self.fig_rzview.axis)
            self.gui.ascot.plotseparatrix(out["r"], out["phi"], out["z"],
                                          out["time"], self.fig_rzview.axis)
            self.fig_rzview.draw()

        self.savebutton.configure(command=storesettings)
        self.plotbutton.configure(command=plot)

        self.panel.pack()
        self.settingsframe.pack()
        self.plotcanvas.pack(fill="both", expand=True)
        plot()
