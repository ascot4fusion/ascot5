"""
Contains definition of OptionsFrame class.

File: optionsframe.py
"""
import tkinter

import numpy as np

class OptionsFrame(tkinter.Frame):
    """
    Frame for viewing simulation options.
    """

    def __init__(self, gui, options):
        super().__init__(gui._root)
        self._gui = gui

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
        self.options   = options
        self.viewfull()


    def _backtoindex(self):
        from .indexframe import IndexFrame
        self._gui.displayframe(IndexFrame(self._gui))


    def viewfull(self):
        opt, info = self.options.read(info=True)

        def trim(val):
            if np.abs(val) >= 1e4 or (np.abs(val) <=1e-4 and val != 0):
                return "{0:e}".format(val)
            else:
                return "{0:g}".format(val)

        t = ""
        for i in info:
            if len(i) == 1:
                t = "\n" + t + i[0] + "\n\n"
                continue

            k = i[0]
            d = i[1]
            v = opt[k]
            t = t + d
            if(len(v) == 1):
                t = t + k + "=" + trim(v[0]) + "\n\n"
            elif(len(v) > 1):
                t = t + k + "=" \
                    + ",".join([trim(vv) for vv in v]) \
                    + "\n\n"

        # Remove empty lines from beginning and end
        t = t[5:-2]

        self.txt.delete("1.0", tkinter.END)
        self.txt.insert("end", t)
        self.txt.configure(state="disabled")
