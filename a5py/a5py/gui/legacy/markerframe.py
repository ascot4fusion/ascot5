"""
Contains definition of MarkerFrame class.

File: markerframe.py
"""
import tkinter

import numpy as np

class MarkerFrame(tkinter.Frame):
    """
    Frame for viewing simulation options.
    """

    def __init__(self, gui, marker):
        super().__init__(gui._root)
        self._gui = gui
        self.marker = marker

        sidepanel  = tkinter.Frame(self, height=gui.height, width=250)
        backbutton = tkinter.Button(sidepanel, text="Back")

        backbutton.config(command=self._backtoindex)
        backbutton.pack()

        sidepanel.grid(row=0, column=0, sticky="NEWS")
        self.grid_propagate(0)

        self.txt = tkinter.Text(self, borderwidth=3, relief="sunken")
        self.txt.grid(row=0, column=1, sticky="nsew", padx=2, pady=2)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        scrollb = tkinter.Scrollbar(self, command=self.txt.yview)
        scrollb.grid(row=0, column=2, sticky='nsew')
        self.txt['yscrollcommand'] = scrollb.set


        self.sidepanel = sidepanel
        self.viewfull()


    def _backtoindex(self):
        from .indexframe import IndexFrame
        self._gui.displayframe(IndexFrame(self._gui))


    def viewfull(self):
        data = self.marker.read()

        text = ""
        n = int(data["n"])

        text = "Number of markers: " + str(n) + "\n\n"
        del data["n"]
        keys = list(data.keys())
        keys.sort()
        for k in keys:
            text += "%12s" % k

        text += "\n"
        for i in range(n):
            for k in keys:
                text += "   %+1.2E" % data[k][i]
            text += "\n"

        self.txt.delete("1.0", tkinter.END)
        self.txt.insert("end", text)
        self.txt.configure(state="disabled")
