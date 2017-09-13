import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import tkinter as tk
from tkinter import ttk

import a5py.ascot5io.ascot5 as ascot5
import a5py.postprocessing.state as stateplot

LARGE_FONT= ("Verdana", 12)


class ascot5GUI(tk.Tk):

    def __init__(self, a5fn, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "ASCOT5 GUI")
        
        self.a5fn = a5fn

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (FrameMain, FrameInputs, FrameInistate, FrameEndstate):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(FrameMain)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()

        
class FrameMain(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Main", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = ttk.Button(self, text="Inputs",
                            command=lambda: controller.show_frame(FrameInputs))
        button.pack()

        button2 = ttk.Button(self, text="Inistate",
                             command=lambda: controller.show_frame(FrameInistate))
        button2.pack()

        button3 = ttk.Button(self, text="Endstate",
                             command=lambda: controller.show_frame(FrameEndstate))
        button3.pack()


class FrameInputs(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Inputs", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Main",
                            command=lambda: controller.show_frame(FrameMain))
        button1.pack()


class FrameInistate(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Inistate", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Main",
                            command=lambda: controller.show_frame(FrameMain))
        button1.pack()

        s = ascot5.read_hdf5(controller.a5fn, "states")["states"]["inistate"]
        stateplot.gather(s)

        f = Figure(figsize=(5,5), dpi=100)
        a = f.add_subplot(2,3,1)
        stateplot.plotEkin(s["Ekin"],ax=a)
        a = f.add_subplot(2,3,2)
        stateplot.plotpitch(s["pitch"],ax=a)
        a = f.add_subplot(2,3,3)
        stateplot.plottime(s["time"],ax=a)
        a = f.add_subplot(2,3,4)
        stateplot.plotrho(s["rho"],ax=a)
        a = f.add_subplot(2,3,5)
        stateplot.plotendcond(s["endCond"],ax=a)
        a = f.add_subplot(2,3,6)
        stateplot.plotRz(s["R"],s["z"],ax=a)

        canvas = FigureCanvasTkAgg(f, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class FrameEndstate(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Endstate", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Main",
                             command=lambda: controller.show_frame(FrameMain))
        button1.pack()

        s = ascot5.read_hdf5(controller.a5fn, "states")["states"]["endstate"]
        stateplot.gather(s)

        f = Figure(figsize=(5,5), dpi=100)
        a = f.add_subplot(2,3,1)
        stateplot.plotEkin(s["Ekin"],ax=a)
        a = f.add_subplot(2,3,2)
        stateplot.plotpitch(s["pitch"],ax=a)
        a = f.add_subplot(2,3,3)
        stateplot.plottime(s["time"],ax=a)
        a = f.add_subplot(2,3,4)
        stateplot.plotrho(s["rho"],ax=a)
        a = f.add_subplot(2,3,5)
        stateplot.plotendcond(s["endCond"],ax=a)
        a = f.add_subplot(2,3,6)
        stateplot.plotRz(s["R"],s["z"],ax=a)

        canvas = FigureCanvasTkAgg(f, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def run(a5fn):
    app = ascot5GUI(a5fn)
    app.mainloop()
