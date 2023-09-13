import tkinter as tk
from tkinter import ttk

from .interactiveorbit import OrbitFrame, OrbitCanvas

class ContentInteractive:

    def __init__(self, gui, settings, canvas):
        self.gui    = gui
        self.canvas = canvas

        runselection = ttk.Notebook(settings)

        runpoincare = ttk.Frame(runselection)
        runorbit    = OrbitFrame(runselection)

        #runselection.add(runpoincare, text="Poincaré")
        runselection.add(runorbit,    text="Orbit")

        canvaspoincare = ttk.Frame(canvas)
        canvasorbit    = OrbitCanvas(canvas)

        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            self.view(tab)

        runselection.bind('<<NotebookTabChanged>>', on_tab_change)
        runselection.pack(fill="both", expand=True)

        runorbit.init(canvasorbit)
        canvasorbit.init()

        self.runpoincare    = runpoincare
        self.runorbit       = runorbit
        self.canvaspoincare = canvaspoincare
        self.canvasorbit    = canvasorbit

        self.active_canvas = canvasorbit

    def display(self):
        pass
        #self.pack()

    def view(self, tab):
        if tab == "Orbit":
            self.gui.ascot.input_free()
            self.active_canvas.pack_forget()
            self.runorbit.display(self.gui)
            self.active_canvas = self.runorbit.canvas
            self.active_canvas.pack(fill="both", expand=True)
