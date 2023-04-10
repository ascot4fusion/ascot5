import tkinter as tk
from tkinter import ttk
import numpy as np

from .components import PlotFrame, NumEntry
from a5py.ascotpy.libbfield  import LibBfield
from a5py.ascotpy.libefield  import LibEfield
from a5py.ascotpy.libplasma  import LibPlasma
from a5py.ascotpy.libneutral import LibNeutral
from a5py.ascotpy.libboozer  import LibBoozer
from a5py.ascotpy.libmhd     import LibMhd

class ContentInput:

    def display_input(self, gui, panel, canvas):

        ## Ascot stuff
        quantities = []

        bfield = False
        if hasattr(gui.ascot.hdf5, "bfield"):
            quantities += LibBfield.quantities
            bfield = True
        efield = False
        if hasattr(gui.ascot.hdf5, "efield"):
            quantities += LibEfield.quantities
            efield = True
        plasma = False
        if hasattr(gui.ascot.hdf5, "plasma"):
            quantities += LibPlasma.quantities
            plasma = True
        neutral = False
        if hasattr(gui.ascot.hdf5, "neutral"):
            quantities += LibNeutral.quantities
            neutral = True
        boozer = False
        if hasattr(gui.ascot.hdf5, "boozer"):
            quantities += LibBoozer.quantities
            boozer = True
        mhd = False
        if hasattr(gui.ascot.hdf5, "mhd"):
            quantities += LibMhd.quantities
            mhd = True

        self.gui.ascot.init(bfield=bfield, efield=efield, neutral=neutral,
                            plasma=plasma, boozer=boozer, mhd=mhd,
                            ignorewarnings=True)

        ## First add widgets ##

        fig_rzview = PlotFrame(canvas)
        fig_rzview.place(relheight=0.8, anchor="nw")

        xmin_entry = NumEntry(panel, labeltext="R = ", defval=0.1)
        xmax_entry = NumEntry(panel, labeltext=" - ",  defval=10)
        xnum_entry = NumEntry(panel, labeltext=" : ",  defval=50,
                                   isint=True)

        ymin_entry = NumEntry(panel, labeltext="z = ", defval=-8)
        ymax_entry = NumEntry(panel, labeltext=" - ",  defval=8)
        ynum_entry = NumEntry(panel, labeltext=" : ",  defval=50,
                                   isint=True)

        xmin_entry.grid(row=0, column=0, sticky="W")
        xmax_entry.grid(row=0, column=1, sticky="W")
        xnum_entry.grid(row=0, column=2, sticky="W")

        ymin_entry.grid(row=1, column=0, sticky="W")
        ymax_entry.grid(row=1, column=1, sticky="W")
        ynum_entry.grid(row=1, column=2, sticky="W")

        phi_entry = NumEntry(panel, labeltext="phi [deg]:", defval=0)
        phi_entry.grid(row=2, column=0)

        time_entry = NumEntry(panel, labeltext="time [s]:", defval=0)
        time_entry.grid(row=2, column=1)

        cmin_entry = NumEntry(panel, labeltext="Color range", defval="")
        cmax_entry = NumEntry(panel, labeltext=" - ", defval="")
        cmin_entry.grid(row=3, column=0)
        cmax_entry.grid(row=3, column=1)

        qchoice = tk.StringVar(panel)
        qchoice.set(quantities[0])
        qinput = ttk.Combobox(panel, width=10, textvariable=qchoice)
        qinput["values"] = quantities
        qinput.grid(row=4, column=1, sticky="WE")

        def plot():

            r = np.linspace( xmin_entry.getval() ,
                             xmax_entry.getval() ,
                             xnum_entry.getval() )

            z = np.linspace( ymin_entry.getval() ,
                             ymax_entry.getval() ,
                             ynum_entry.getval() )

            phi  = phi_entry.getval() * np.pi / 180
            time = time_entry.getval()

            clim = [None, None]
            if not cmin_entry.isempty():
                clim[0] = cmin_entry.getval()
            if not cmax_entry.isempty():
                clim[1] = cmax_entry.getval()

            fig_rzview.clear()
            gui.ascot.plotRz(r, phi, z, time, qchoice.get(),
                             axes=fig_rzview.axis, clim=clim)
            gui.ascot.plotseparatrix(r, phi, z, time, fig_rzview.axis)
            fig_rzview.draw()


        plotbutton = tk.Button(panel, text="Plot", command=plot)
        plotbutton.grid(row=4, column=2, sticky="WE")
