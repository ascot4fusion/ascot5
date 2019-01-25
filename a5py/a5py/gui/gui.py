"""
GUI main module.

File gui.py
"""

import os
import sys
import tkinter
from tkinter.filedialog import askopenfilename, askdirectory
from tkinter import messagebox
import a5py.ascot5io.ascot5 as ascot5

from .indexframe import IndexFrame

class GUI:
    """
    GUI object that controls the window and has direct access to the Ascot data.

    The ascot file is kept constantly open, so whenever the data in HDF5 file
    has changed, GUI.reload() must be called which reads the file and
    re-initializes GUI.
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
            self.ask_openfile()

        else:
            self._h5fn  = os.path.abspath(h5fn)
            self._ascot = ascot5.Ascot(self._h5fn)

        self._root.protocol("WM_DELETE_WINDOW", self.close)

        self._ascotfolder = ""
        self._ascotpy     = None

    def get_ascotobject(self):
        """
        Get Ascot object of the currently opened HDF5 file.
        """
        return self._ascot

    def get_ascotfilename(self):
        """
        Get filename of currently opened HDF5 file.
        """
        return self._h5fn

    def get_ascotfolder(self):
        """
        Get filename of currently opened HDF5 file.
        """
        return self._ascotfolder

    def get_ascotpy(self):
        """
        Get interface to a running Ascot process.
        """
        return self._ascotpy

    def get_root(self):
        """
        Get root windonw for displaying frames.
        """
        return self._root

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

        Args:
            newframe : Frame <br>
                The frame to be displayed.
        """
        if self._current is not None:
            self._current.pack_forget()
            self._current.destroy()

        self._current = newframe
        self._current.pack_propagate(0)
        self._current.pack(fill=tkinter.BOTH, expand=1)

    def ask_openfile(self):
        """
        Open dialog for opening filename and open it.
        """
        fn = askopenfilename( title="Select file",
                              filetypes = [("HDF5 files","*.h5")] )
        if len(fn) == 0:
            pass
        else:
            self._h5fn  = os.path.abspath(fn)
            self._ascot = ascot5.Ascot(self._h5fn)

    def ask_openfolder(self):
        """
        Open dialog for opening a folder where Ascot source code is.
        """
        fn = askdirectory(title="Select folder")

        if (len(fn) == 0) :
            pass

        elif not os.path.isfile(fn + "/ascotpy.so") or \
             not os.path.isfile(fn + "/ascotpymod.py"):
            messagebox.showerror("Failed to open folder",
                                 "Folder does not contain ascotpy.so "
                                 + "or ascotpymod.py")
        else:
            fn = os.path.abspath(fn)
            sys.path.insert(0, fn)
            from ascotpymod import Ascotpy
            self._ascotfolder = fn
            self._ascotpy = Ascotpy(fn + "/ascotpy.so", self._h5fn)

    def reload(self):
        """
        Reinitializes GUI (if Ascot data has changed or different file opened).
        """
        self._ascot = ascot5.Ascot(self._h5fn)
        if self._ascotpy is not None:
            self._ascotpy.ascotpy_reload(self._h5fn)
        self.displayframe( IndexFrame(self) )

    def close(self):
        """
        Close the gui and terminate the program.
        """
        self._root.destroy()
        exit()
