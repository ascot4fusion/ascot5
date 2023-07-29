import tkinter as tk
from tkinter import ttk

import numpy as np
import unyt

from .components import ContentTab, PlotFrame, NumEntry

class Field(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)

        plotrz = PlotFrame(canvas)
        self.canvas.slideadd(self, plotrz)

        f1 = ttk.Frame(self)
        f2 = ttk.Frame(self)

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

        def storesettings():
            gui.params.store(
                gui.ascot.data,
                ["input_rzplot_minr", "input_rzplot_maxr", "input_rzplot_numr",
                 "input_rzplot_minz", "input_rzplot_maxz", "input_rzplot_numz",
                 "input_rzplot_qnt"] )

        def plot():
            r = np.linspace( xmin_entry.getval() ,
                             xmax_entry.getval() ,
                             xnum_entry.getval() )

            z = np.linspace( ymin_entry.getval() ,
                             ymax_entry.getval() ,
                             ynum_entry.getval() )

            r *= unyt.m
            z *= unyt.m
            phi  = phi_entry.getval()  * unyt.deg
            time = time_entry.getval() * unyt.s

            clim = [None, None]
            if not cmin_entry.isempty():
                clim[0] = cmin_entry.getval()
            if not cmax_entry.isempty():
                clim[1] = cmax_entry.getval()
                qnt  = qchoice.get()

            qnt = self.qchoice.get()

            plotrz.clear()
            self.gui.ascot.input_plotrz(r, z, qnt, phi=phi, t=time, clim=clim,
                                        axes=plotrz.axes)
            #self.gui.ascot.input_plotseparatrix(r, phi, z, time,
            #                                    axes=plorrz.axes)
            plotrz.draw()

        savebutton.configure(command=storesettings)
        plotbutton.configure(command=plot)

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

    def selecttab(self):
        qnt, desc = zip(*self.gui.ascot.input_eval_list(show=False).items())
        self.qinput["values"] = qnt
        if self.qchoice.get() is None or self.qchoice.get() not in qnt:
            self.qchoice.set(qnt[0])

        self.canvas.slideshow(self)

class Preflight(ContentTab):

    def __init__(self, frame, canvas, gui):
        super().__init__(frame, canvas, gui)
        self.text = tk.Text(self, height=25, width=50)
        self.text.pack()
        self.plot = Preflight.Canvas(canvas)
        self.canvas.slideadd(self, self.plot)

    def selecttab(self):
        """Run and display results of the preflight tests.
        """
        msg = []
        msg += self.gui.ascot.preflight_inputspresent()
        msg += self.gui.ascot.preflight_optionsconsistent()
        #msg += check_bfield_psi0(gui.ascot)
        if len(msg) == 0:
            msg += ["Preflight checks completed with no issues."]

        msg0 = ""
        for m in msg:
            msg0 += m + "\n"

        self.text.configure(state='normal')
        self.text.delete("1.0", "end")
        self.text.insert("end", msg0)
        self.text.configure(state='disabled')

        self.plot.clear()
        self.gui.ascot.preflight_plottopview(axes=self.plot.axes[0])
        #plot_energypitch(ascot, axes=self.axes[1])
        #plot_rhophi(ascot, axes=self.axes[2])
        self.plot.draw()
        self.canvas.slideshow(self)

    class Canvas(PlotFrame):

        def set_axes(self):
            ax1 = self.fig.add_subplot(2,2,(1,3))
            ax2 = self.fig.add_subplot(2,2,2)
            ax3 = self.fig.add_subplot(2,2,4)
            return [ax1, ax2, ax3]
