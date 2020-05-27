import os

import tkinter as tk
from tkinter import ttk

from .optionsmanager import OptionsManager

class ContentManager(OptionsManager):
    """
    Manages the contents in SettingsFrame and in CanvasFrame.
    """


    def __init__(self, settingsframe, canvasframe):

        # These two frames are kept immutable
        self.settingsframe_original = settingsframe
        self.canvasframe_original   = canvasframe

        # These two form an extra layer and they can be changed
        self.settingsframe = tk.Frame(self.settingsframe_original)
        self.canvasframe   = tk.Frame(self.canvasframe_original)

        # (The frames are initialized here)
        self.clear()

        # Show ASCOT5 logo at startup.
        canvas = tk.Canvas(self.settingsframe)
        canvas.pack(expand=True, fill="both")

        # self.logo prevents image being cleared by garbage collector
        logo = os.path.join(os.path.dirname(__file__), "logo.png")
        self.logo = tk.PhotoImage(file=logo)
        canvas.create_image(150, 50, image=self.logo)


    def clear(self):
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


    def selectionchanged(self, parent, qid, ascot, ascotpy):
        self.clear()
        if qid is None:
            # For parent groups there is nothing to display.
            return

        if parent == "results":
            group = ascot["q"+qid]
        else:
            group = ascot[parent]["q"+qid]

        # Always show description box on top.
        f1 = tk.Frame(self.settingsframe)
        tk.Label(f1, text="Description:").pack(side="left")
        f1.grid(row=0, column=0, sticky="ew")

        f2 = tk.Frame(self.settingsframe)
        save = tk.Button(f2, text="Save")
        save.pack(side="right")

        load = tk.Button(f2, text="Revert")
        load.pack(side="right")

        f2.grid(row=0, column=1, sticky="ew")

        descbox = tk.Text(self.settingsframe, height=5, width=50)
        descbox.grid(row=1, column=0, columnspan=2, sticky="nsew",
                     padx=2, pady=2)
        descbox.insert("end", group.get_desc())

        # Helper function
        def replacetext(text):
            descbox.delete("1.0", "end")
            descbox.insert("end", text)

        # Button functionality
        save.configure(
            command=lambda:group.set_desc(descbox.get("1.0", "end-1c")))
        load.configure(
            command=lambda:replacetext(group.get_desc()))

        # Add extra frame where widgets can be added
        frame = tk.Frame(self.settingsframe)
        frame.grid(row=2, column=0, columnspan=2, sticky="nsew")

        if parent == "options":
            self.activateoptions(frame, self.canvasframe, group, ascotpy)
