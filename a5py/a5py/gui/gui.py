"""
GUI main module.

File core.py
"""

import os
import tkinter
from tkinter.filedialog import askopenfilename
import a5py.ascot5io.ascot5 as ascot5

from .indexframe import IndexFrame

class GUI:
    """
    GUI object that controls the window and has direct access to the Ascot data.
    """

    def __init__(self, h5fn=None):
        """
        Spawn a new GUI.

        During initialization, window is initialized (but not shown) and Ascot
        object is initalized from the filename. If filename was not given, a
        window prompting the user to choose a file is shown.
        """
        self._root = tkinter.Tk()
        self._root.withdraw()

        self._current = None

        sw = self._root.winfo_screenwidth()
        sh = self._root.winfo_screenheight()

        # This is ugly, but makes window work if there is a dual monitor
        if sw > 2000:
            sw = 2000

        w = sw*3/4
        h = sh*3/4
        x = (sw/2) - (w/2)
        y = (sh/2) - (h/2)

        self.width  = w
        self.height = h
        self._root.geometry('%dx%d+%d+%d' % (w, h, x, y))
        self._root.resizable(0, 0)
        self._root.title("ASCOT5 GUI")

        if h5fn is None:
            h5fn = askopenfilename(title="Choose a file to open.")
            self._h5fn  = os.path.abspath(h5fn)
            self._ascot = ascot5.Ascot(self._h5fn)

        else:
            self._h5fn  = os.path.abspath(h5fn)
            self._ascot = ascot5.Ascot(self._h5fn)

    def launch(self):
        """
        Show GUI and display index frame.
        """
        self._root.deiconify()
        self.displayframe( IndexFrame(self) )
        self._root.mainloop()

    def displayframe(self, newframe):
        """
        Display a new frame and delete the old one.
        """
        if self._current is not None:
            self._current.pack_forget()
            self._current.destroy()

        self._current = newframe
        self._current.pack_propagate(0)
        self._current.pack(fill=tkinter.BOTH, expand=1)

    def reload(self):
        """
        Reinitializes GUI (if Ascot data has changed or different file opened).
        """
        self._ascot = ascot5.Ascot(self._h5fn)
        self.displayframe( IndexFrame(self) )


if __name__ == '__main__':
    # Make this script callable for development purposes.
    import sys
    gui = GUI(sys.argv[1])
    gui.launch()
