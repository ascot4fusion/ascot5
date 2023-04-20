import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class SummaryFrame(ttk.Frame):

    def init(self, canvas):
        self.canvas = canvas

        return self


    def display(self):
        pass


class SummaryCanvas(ttk.Frame):

    def init(self):
        return self
