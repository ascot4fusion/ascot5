import os

import tkinter as tk
from tkinter import ttk

from .contentgroup       import ContentGroup
from .contentprecheck    import ContentPrecheck
from .contentinput       import ContentInput
from .contentoutput      import ContentOutput
from .contentinteractive import ContentInteractive

class ContentManager(ContentGroup, ContentPrecheck, ContentInput,
                     ContentOutput, ContentInteractive):
    """
    Manages the contents in SettingsFrame and in CanvasFrame.
    """


    def __init__(self, gui, settingsframe, canvasframe):

        self.gui = gui

        # These two frames are kept immutable
        self.settingsframe_original = settingsframe
        self.canvasframe_original   = canvasframe

        # These two form an extra layer and they can be changed
        self.settingsframe = tk.Frame(self.settingsframe_original)
        self.canvasframe   = tk.Frame(self.canvasframe_original)

        # Load graphics
        logo = os.path.join(os.path.dirname(__file__), "logo.png")
        self.logo = tk.PhotoImage(file=logo)

        # (The frames are initialized here)
        self.clear_content()


    def clear_content(self):
        """
        Clear settings and canvas frames.
        """
        self.settingsframe.pack_forget()
        self.canvasframe.pack_forget()
        self.settingsframe.destroy()
        self.canvasframe.destroy()

        self.settingsframe = tk.Frame(self.settingsframe_original)
        self.settingsframe.pack(fill="both", expand=True)

        self.canvasframe = tk.Frame(self.canvasframe_original, bg="white")
        self.canvasframe.pack(fill="both", expand=True)


    def display_content(self, content, **kwargs):
        """
        Parse what content to display and display it.

        All display functions should be accessed via this interface.
        """

        if content == "welcome":
            self.clear_content()
            self._display_welcome(
                self.gui, self.settingsframe, self.canvasframe
            )

        if content == "group":
            self.clear_content()
            self._display_group(
                self.gui, self.settingsframe, self.canvasframe,
                kwargs["parent"], kwargs["qid"]
            )

        if content == "preflight":
            self.clear_content()
            self.gui.ascot.init(bfield=True, ignorewarnings=True)
            self.display_precheck(
                self.gui, self.settingsframe, self.canvasframe
            )
        if content == "input":
            self.clear_content()
            self.display_input(
                self.gui, self.settingsframe, self.canvasframe
            )
        if content == "output":
            self.clear_content()
            self.gui.ascot.init_from_run(self.gui.ascot.hdf5.active,
                bfield=True, efield=True, neutral=True,
                plasma=True, boozer=True, mhd=True)
            self.display_output(
                self.gui, self.settingsframe, self.canvasframe,
                self.gui.ascot.hdf5.active
            )
        if content == "interactive":
            pass


    def _display_welcome(self, gui, settingsframe, canvasframe):
        """
        Display welcome screen.
        """

        # Show ASCOT5 logo
        canvas = tk.Canvas(self.settingsframe)
        canvas.pack(expand=True, fill="both")
        canvas.create_image(150, 50, image=self.logo)
