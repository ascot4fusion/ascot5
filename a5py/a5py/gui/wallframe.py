"""
Contains definition of WallFrame class.

File: wallframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from .plotframe import PlotFrame
from .components import NumEntry
from .cameraControlPanel import generateCameraControlPanel

class WallFrame(PlotFrame):
    """
    A frame for plotting wall inputs.
    """

    def __init__(self, gui, wall):
        """
        Initialize and show default plot.
        """

        super().__init__(gui)
        self._gui  = gui
        self._wall = wall

        panel = self.get_sidepanel()
        self.phi_entry = NumEntry(panel, labeltext="phi [deg]:", defval=0)
        self.phi_entry.pack() #.grid(row=0, column=0)

        plotbutton = tkinter.Button(panel, text="Plot", command=self._plot)
        plotbutton.pack() #.grid(row=1, column=0, sticky="WE")


        vtkText=tkinter.Label(panel, text="Plot the 3D wall in separate window:")
        vtkText.pack() #.grid(row=2, column=0, sticky="WE")
        
        vtkbutton = tkinter.Button(panel, text="Plot in 3D with VTK", command=self._plotVTKwall)
        
        self.camControlPanelComponents,self.applyCameraControlPanel= generateCameraControlPanel(panel)
        self.camControlPanelComponents['panel'].pack()

        
        vtkbutton.pack() #.grid(row=3, column=0, sticky="WE")



        self._plot()


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        phi  = self.phi_entry.getval() * np.pi / 180

        fig  = self.get_fig()
        axes = fig.add_subplot(1,1,1)
        self._wall.plotRz(axes=axes, phi=phi)

        self.draw()

    def _plotVTKwall(self):
        
        import a5py.wall.a5vtkwall
       
        if not hasattr(self, 'vtkWall'):
            print("Starting to generate VTK 3D-wall...") 
            self.vtkWall=a5py.wall.a5vtkwall.a5VtkWall(W=self._wall)
            print("Add triangle numbering for colors")
            self.vtkWall.addIndex()
            print(" ...done. Next plotting") 
        
        self.camControl = a5py.wall.a5vtkwall.camControl()
        
        self.applyCameraControlPanel(self.camControlPanelComponents,self.camControl)

        print('Using camControl: '+str(self.camControl))

        self.vtkWall.plot(logarithmicColorScale=False,camControl=self.camControl)
       

