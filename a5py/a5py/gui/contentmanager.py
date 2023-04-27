"""
Defines ContentManager class for showing stuff at Canvas and Settings frames.

File: contentmanager.py
"""
import tkinter as tk
from tkinter import ttk

from .contentgroup       import ContentGroup
from .contentprecheck    import ContentPrecheck
from .contentinput       import ContentInput
from .contentoutput      import ContentOutput
from .contentinteractive import ContentInteractive

class ContentManager():
    """
    Manages the contents in SettingsFrame and in CanvasFrame.

    Settings frame is a notebook widget where changing the tab also changes
    what is shown on Canvas*. This is done so that we always pack_forget the
    current frame on Canvas and pack the new one (Note that the frames are
    kept as that way we preserve the settings user has given previously, which
    makes smoother user-experience when different frames are toggled).

    *For group contents the contents on Canvas also depends what is shown in
    the group treeview.

    This is fairly simple class as the actual contents are defined in Content*
    classes that are imported an used here. If those also have subcontents and
    implement a notebook to choose what is display, then those classes mimic the
    structure in this class. In theory, everything could be defined here but
    that would make this class huge...
    """


    def __init__(self, gui, settingsframe, canvasframe):
        """
        Generate all content widgets (which also generate their contents).
        """

        self.gui = gui

        # Add tabs to the notebook widget
        settings         = ttk.Notebook(settingsframe)
        groupframe       = ttk.Frame(settings)
        preflightframe   = ttk.Frame(settings)
        inputframe       = ttk.Frame(settings)
        outputframe      = ttk.Frame(settings)
        interactiveframe = ttk.Frame(settings)
        settings.add(groupframe,       text="Group")
        settings.add(preflightframe,   text="Preflight")
        settings.add(inputframe,       text="Input")
        settings.add(outputframe,      text="Analysis")
        settings.add(interactiveframe, text="Run")

        # Initialize content frames (nothing is shown yet)
        groupcanvas       = ttk.Frame(canvasframe)
        preflightcanvas   = ttk.Frame(canvasframe)
        inputcanvas       = ttk.Frame(canvasframe)
        outputcanvas      = ttk.Frame(canvasframe)
        interactivecanvas = ttk.Frame(canvasframe)
        self.contentgroup        = ContentGroup(
            gui, groupframe, groupcanvas)
        self.contentprecheck    = ContentPrecheck(
            gui, preflightframe, preflightcanvas)
        self.contentinput       = ContentInput(
            gui, inputframe, inputcanvas)
        self.contentoutput      = ContentOutput(
            gui, outputframe, outputcanvas)
        self.contentinteractive = ContentInteractive(
            gui, interactiveframe, interactivecanvas)

        # Have an empty canvas initially
        self.active_canvas = ttk.Frame(canvasframe)
        self.active_canvas.pack(fill="both", expand=True)

        # Display contents when tab changes
        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            self.display_content(tab)

        settings.bind('<<NotebookTabChanged>>', on_tab_change)
        settings.pack(fill="both", expand=True)


    def update_content(self):
        """
        Redraw contents on settings and canvas frames.
        """
        self.display_content(self.content)


    def display_content(self, content):
        """
        Parse what content to display and display it.

        All display functions should be accessed via this interface.
        """
        self.content = content

        if content == "Group":
            self.active_canvas.pack_forget()
            self.contentgroup.display()
            self.active_canvas = self.contentgroup.canvas
            self.active_canvas.pack(fill="both", expand=True)

        if content == "Preflight":
            self.active_canvas.pack_forget()
            self.contentprecheck.display()
            self.active_canvas = self.contentprecheck.canvas
            self.active_canvas.pack(fill="both", expand=True)

        if content == "Input":
            self.active_canvas.pack_forget()
            self.contentinput.display()
            self.active_canvas = self.contentinput.canvas
            self.active_canvas.pack(fill="both", expand=True)

        if content == "Analysis":
            self.active_canvas.pack_forget()
            self.contentoutput.display()
            self.active_canvas = self.contentoutput.canvas
            self.active_canvas.pack(fill="both", expand=True)

        if content == "Run":
            self.active_canvas.pack_forget()
            self.contentinteractive.display()
            self.active_canvas = self.contentinteractive.canvas
            self.active_canvas.pack(fill="both", expand=True)
