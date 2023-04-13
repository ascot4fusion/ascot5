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
    """


    def __init__(self, gui, settingsframe, canvasframe):

        self.gui = gui

        # Add tabs to the notebook widget
        groupframe       = ttk.Frame(canvasframe)
        preflightframe   = ttk.Frame(canvasframe)
        inputframe       = ttk.Frame(canvasframe)
        outputframe      = ttk.Frame(canvasframe)
        interactiveframe = ttk.Frame(canvasframe)
        canvasframe.add(groupframe,       text="Group")
        canvasframe.add(preflightframe,   text="Preflight")
        canvasframe.add(inputframe,       text="Input")
        canvasframe.add(outputframe,      text="Analysis")
        canvasframe.add(interactiveframe, text="Run")

        # Initialize content frames (nothing is shown yet)
        self.contentgroup       = ContentGroup(
            gui, settingsframe, groupframe)
        self.contentprecheck    = ContentPrecheck(
            gui, settingsframe, preflightframe)
        self.contentinput       = ContentInput(
            gui, settingsframe, inputframe)
        self.contentoutput      = ContentOutput(
            gui, settingsframe, outputframe)
        self.contentinteractive = ContentInteractive(
            gui, settingsframe, interactiveframe)

        # Have an empty frame here initially
        self.active_settingsframe = ttk.Frame(settingsframe)
        self.active_settingsframe.pack(fill="both")

        # Display contents when tab changes
        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            self.display_content(tab)

        canvasframe.bind('<<NotebookTabChanged>>', on_tab_change)


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
            self.active_settingsframe.pack_forget()
            tree = self.gui.groups.tree
            qid    = tree.item(tree.selection(), "text")
            parent = tree.item(tree.parent(tree.selection()), "text")
            self.contentgroup.display(parent, qid)
            self.active_settingsframe = self.contentgroup.settingsframe

        if content == "Preflight":
            self.active_settingsframe.pack_forget()
            self.contentprecheck.display()
            self.active_settingsframe = self.contentprecheck.settingsframe

        if content == "Input":
            self.active_settingsframe.pack_forget()
            self.contentinput.display()
            self.active_settingsframe = self.contentinput.settingsframe

        if content == "Analysis":
            self.active_settingsframe.pack_forget()
            self.contentoutput.display()
            self.active_settingsframe = self.contentoutput.settingsframe

        if content == "Run":
            self.active_settingsframe.pack_forget()
            self.contentinteractive.display()
            self.active_settingsframe = self.contentinteractive.settingsframe
