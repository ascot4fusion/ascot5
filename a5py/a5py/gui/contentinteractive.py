import tkinter as tk
from tkinter import ttk

class ContentInteractive:

    def __init__(self, gui, settings, canvas):
        self.gui    = gui
        self.canvas = canvas

        runselection = ttk.Notebook(settings)

        runpoincare = ttk.Frame(runselection)
        runorbit    = ttk.Frame(runselection)

        runselection.add(runpoincare, text="Poincaré")
        runselection.add(runorbit,    text="Orbit")

        canvaspoincare = ttk.Frame(canvas)
        canvasorbit    = ttk.Frame(canvas)

        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            #self.display_content(tab)

        runselection.bind('<<NotebookTabChanged>>', on_tab_change)
        runselection.pack(fill="both", expand=True)

    def display(self):
        pass
