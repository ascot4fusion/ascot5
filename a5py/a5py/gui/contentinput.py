import tkinter as tk
from tkinter import ttk
import numpy as np

from .components import PlotFrame, NumEntry

class ContentInput:

    def __init__(self, gui, settings, canvas):
        self.gui = gui

        ## Add widgets

        class EntryFrame(ttk.Frame):

            def init(self, canvas):
                f = ttk.Frame(self)
                f.pack(fill="both", expand=True)
                f1 = ttk.Frame(f)
                f2 = ttk.Frame(f)

                f1.pack(anchor="nw", side="left")
                f2.pack(anchor="ne", side="left", fill="x", expand=True)


                xmin_entry = NumEntry(f1, labeltext="R [m]     = ", entrywidth=5,
                                      labelwidth=12, anchor="w", defval=0.1)
                xmax_entry = NumEntry(f1, labeltext="–",    entrywidth=5,
                                      labelwidth=2, anchor="c", defval=10.0)
                xnum_entry = NumEntry(f1, labeltext="x",    entrywidth=5,
                                      labelwidth=2, anchor="c", defval=50,
                                      isint=True)

                ymin_entry = NumEntry(f1, labeltext="z [m]     = ", entrywidth=5,
                                      labelwidth=12, anchor="w", defval=8.0)
                ymax_entry = NumEntry(f1, labeltext="–",    entrywidth=5,
                                      labelwidth=2, anchor="c", defval=-8.0)
                ynum_entry = NumEntry(f1, labeltext="x",    entrywidth=5,
                                      labelwidth=2, anchor="c", defval=50,
                                      isint=True)

                xmin_entry.grid(row=0, column=0)
                xmax_entry.grid(row=0, column=1)
                xnum_entry.grid(row=0, column=2)

                ymin_entry.grid(row=1, column=0)
                ymax_entry.grid(row=1, column=1)
                ynum_entry.grid(row=1, column=2)

                phi_entry = NumEntry(f1, labeltext="phi [deg] = ", entrywidth=5,
                                     labelwidth=12, anchor="w", defval=0)
                phi_entry.grid(row=2, column=0)

                time_entry = NumEntry(f1, labeltext="Time [s]  = ", entrywidth=5,
                                      labelwidth=12, anchor="w", defval=0)
                time_entry.grid(row=3, column=0)

                cmin_entry = NumEntry(f1, labeltext="Color     = ", entrywidth=5,
                                      labelwidth=12, anchor="w", defval="")
                cmax_entry = NumEntry(f1, labeltext="-", entrywidth=5,
                                      labelwidth=2, anchor="c", defval="")
                cmin_entry.grid(row=4, column=0)
                cmax_entry.grid(row=4, column=1)

                qchoice = tk.StringVar(f1)
                qinput = ttk.Combobox(f1, width=12, textvariable=qchoice,
                                      state="readonly")
                qinput.grid(row=5,column=0, sticky="W")

                plotbutton = tk.Button(f2, text="Plot", width=3)
                plotbutton.pack(anchor="e")
                savebutton = tk.Button(f2, text="Store", width=3)
                savebutton.pack(anchor="e")

                # Store parameters
                gui.params.add(
                    input_rzplot_minr=xmin_entry.choice,
                    input_rzplot_maxr=xmax_entry.choice,
                    input_rzplot_numr=xnum_entry.choice,
                    input_rzplot_minz=ymin_entry.choice,
                    input_rzplot_maxz=ymax_entry.choice,
                    input_rzplot_numz=ynum_entry.choice,
                    input_rzplot_qnt=qchoice
                )

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
                self.savebutton = savebutton
                self.plotbutton = plotbutton
                self.canvas     = canvas

                return self


            def view(self, gui, quantities):

                self.qinput["values"] = quantities
                if self.qchoice.get() is None or \
                   ( self.qchoice.get() not in quantities ):
                    self.qchoice.set(quantities[0])

                def storesettings():
                    gui.params.store(gui.ascot.h5fn,
                                     ["input_rzplot_minr",
                                      "input_rzplot_maxr",
                                      "input_rzplot_numr",
                                      "input_rzplot_minz",
                                      "input_rzplot_maxz",
                                      "input_rzplot_numz",
                                      "input_rzplot_qnt"]
                                 )
                self.savebutton.configure(command=storesettings)

                def plot():
                    r = np.linspace( self.xmin_entry.getval() ,
                                     self.xmax_entry.getval() ,
                                     self.xnum_entry.getval() )

                    z = np.linspace( self.ymin_entry.getval() ,
                                     self.ymax_entry.getval() ,
                                     self.ynum_entry.getval() )

                    phi  = self.phi_entry.getval() * np.pi / 180
                    time = self.time_entry.getval()

                    clim = [None, None]
                    if not self.cmin_entry.isempty():
                        clim[0] = self.cmin_entry.getval()
                    if not self.cmax_entry.isempty():
                        clim[1] = self.cmax_entry.getval()
                    qnt  = self.qchoice.get()

                    self.canvas.fig_rzview.clear()
                    gui.ascot.plotRz(r, phi, z, time, qnt, clim=clim,
                                     axes=self.canvas.fig_rzview.axes)
                    gui.ascot.plotseparatrix(r, phi, z, time,
                                             self.canvas.fig_rzview.axes)
                    self.canvas.fig_rzview.draw()

                self.plotbutton.configure(command=plot)
                plot()


        class EntryCanvas(ttk.Frame):

            def init(self):
                fig_rzview = PlotFrame(self)
                fig_rzview.place(relheight=0.8, anchor="nw")

                self.fig_rzview = fig_rzview
                return self


        canvasentry = EntryCanvas(canvas).init()
        frameentry = EntryFrame(settings).init(canvasentry)
        canvasentry.pack(fill="both", expand=True)
        frameentry.pack(fill="both", expand=True)

        self.frameentry    = frameentry
        self.canvasentry   = canvasentry
        self.canvas        = canvas


    def display(self):

        # Find what quantities we can plot and init
        quantities = []

        bfield = False
        if hasattr(self.gui.ascot.data, "bfield"):
            quantities += LibBfield.quantities
            bfield = True
        efield = False
        if hasattr(self.gui.ascot.data, "efield"):
            quantities += LibEfield.quantities
            efield = True
        plasma = False
        if hasattr(self.gui.ascot.data, "plasma"):
            quantities += LibPlasma.quantities
            plasma = True
        neutral = False
        if hasattr(self.gui.ascot.data, "neutral"):
            quantities += LibNeutral.quantities
            neutral = True
        boozer = False
        if hasattr(self.gui.ascot.data, "boozer"):
            quantities += LibBoozer.quantities
            boozer = True
        mhd = False
        if hasattr(self.gui.ascot.data, "mhd"):
            quantities += LibMhd.quantities
            mhd = True

        def init():
            self.gui.ascot.init(bfield=bfield, efield=efield, neutral=neutral,
                            plasma=plasma, boozer=boozer, mhd=mhd,
                            ignorewarnings=True)

        msg = self.gui.pleasehold("Ascotpy is being initialized...", init)

        self.frameentry.view(self.gui, quantities)
