"""
Contains definition of DebugFrame class.

File: debugframe.py
"""
import tkinter
import numpy as np

class DebugFrame(tkinter.Frame):
    """
    An debug frame for evaluating data from live action process.
    """

    def __init__(self, gui, ascotpy):
        """
        Initialize index .
        """
        super().__init__(gui._root)
        self._gui     = gui
        self._ascotpy = ascotpy

        self._ascotpy.ascotpy_init(bfield=True)

        indexpanel = tkinter.Frame(self, width=200, height=100)
        toppanel   = tkinter.Frame(self, width=gui.width-100, height=100)
        sidepanel  = tkinter.Frame(self, height=gui.height-100, width=200)
        textpanel  = tkinter.Frame(self)
        backbutton = tkinter.Button(indexpanel, text="Back")

        backbutton.config(command=self._backtoindex)

        backbutton.pack()

        indexpanel.grid(row=0, column=0)
        toppanel.grid(  row=0, column=1)
        sidepanel.grid( row=1, column=0)
        textpanel.grid( row=1, column=1)

        indexpanel.pack_propagate(0)
        toppanel.grid_propagate(0)
        sidepanel.grid_propagate(0)

        self._toppanel  = toppanel
        self._sidepanel = sidepanel

        self._rchoice   = tkinter.StringVar(self)
        self._phichoice = tkinter.StringVar(self)
        self._zchoice   = tkinter.StringVar(self)

        self._rchoice.set(5)
        self._phichoice.set(0)
        self._zchoice.set(0)

        evalbutton = tkinter.Button(sidepanel, text="Evaluate",
                                    command=self._evaluate)

        rlabel   = tkinter.Label(sidepanel, text="R [m]:")
        rentry   = tkinter.Entry(sidepanel, validate = 'key',
                                 width=6, textvariable=self._rchoice)
        philabel = tkinter.Label(sidepanel, text="phi [deg]:")
        phientry = tkinter.Entry(sidepanel, validate = 'key',
                                 width=6, textvariable=self._phichoice)
        zlabel   = tkinter.Label(sidepanel, text="z [m]:")
        zentry   = tkinter.Entry(sidepanel, validate = 'key',
                                 width=6, textvariable=self._zchoice)

        rlabel.grid(    row=0, column=0)
        rentry.grid(    row=0, column=1)
        philabel.grid(  row=1, column=0)
        phientry.grid(  row=1, column=1)
        zlabel.grid(    row=2, column=0)
        zentry.grid(    row=2, column=1)
        evalbutton.grid(row=3, column=1)


        self._output = tkinter.Text(textpanel)
        self._output.pack()
        self._output.configure(state="disabled")



    def _backtoindex(self):
        self._ascotpy.ascotpy_free(bfield=True)
        from .indexframe import IndexFrame
        self._gui.displayframe(IndexFrame(self._gui))


    def _evaluate(self):
        r   = np.array([ float(self._rchoice.get()) ])
        phi = np.array([ float(self._phichoice.get()) * np.pi / 180 ])
        z   = np.array([ float(self._zchoice.get()) ])

        bdata = self._ascotpy.ascotpy_eval_bfield(r, phi, z, evalrho=True,
                                                  evalpsi=True)

        out = ""
        out += "rho : " + str(bdata["rho"]) + "\n"
        out += "psi : " + str(bdata["psi"]) + "\n"

        self._output.configure(state="normal")
        self._output.delete("1.0", tkinter.END)
        self._output.insert("end", out)
        self._output.configure(state="disabled")
