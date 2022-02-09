'''
This module is the interface to VTK for walls.

Created on Mar 6, 2020

@author: sjjamsa
'''

import vtk
import numpy as np



class a5VtkWall(object):
    '''
    This class represents an ASCOT5 wall in VTK polydata datatype.
    The data can include zero or more scalar values per polygon.
    '''


    def __init__(self, W=None):
        '''
        @param W: (optional) The wall_3D class:  a5py.ascot5io.ascot5.Ascot(filename).active.wall
        '''
        
        if W is not None:
            self.fromA5wall(W)
            return
    
        
    def fromA5wall(self, W):
        '''
        @param W: The wall_3D class:  a5py.ascot5io.ascot5.Ascot(filename).active.wall
        '''
        
        points,vertices  = W.getAspointsAndVertices()
        
        self.fromPointsAndVertices(points,vertices)
    
    def fromPointsAndVertices(self,points,vertices):
        '''
        @param points: nx3 array of triangle corners; a list of (x,y,z)-tuples
        @param vertices: nx3 array of 0-based indexes in the points list 
        '''
        
        
        print('Creating VTK points from triangle vertices.')

        
        nTriangles = vertices.shape[0]
        nPoints    = points.shape[0]
        print('Points: {}, Triangles {}'.format(nPoints,nTriangles))
        
        Points = vtk.vtkPoints()
        Points.SetDataTypeToFloat()
        #Points.SetDataTypeToDouble()
        Points.SetNumberOfPoints(nPoints) 
        for i,P in enumerate(points):
            Points.SetPoint(i,P[0], P[1], P[2])



            
        
        print('Creating VTK triangles referencing the points.')
        
        
        Triangles = vtk.vtkCellArray()
        
        # Prepare a connectivity array
        Cells= vtk.vtkIdTypeArray()
        Cells.SetNumberOfTuples(4*nTriangles)
        
        # Fill the array with nVertices, vert1, vert2, ..., vertn, nVertices, vert1, vert2, ...        
        for iTri in range(nTriangles):
            ind=iTri*4
            Cells.SetValue(ind,3) # triangle (3 vertices)
            for iP in range(3):
                ind = iTri*4+iP+1
                Cells.SetValue(ind,vertices[iTri,iP])


        Triangles.SetCells(   nTriangles,   Cells ) 

        
        print('Combining the points and triangles to a single VTK PolyData datastructure.')
        
            
        Polydata = vtk.vtkPolyData()
        Polydata.SetPoints(Points)
        Polydata.SetPolys(Triangles)
        
        Polydata.Modified()
        
        self.nTriangles = nTriangles
        self.vtkPolyData = Polydata

    def addFlag(self,flag,flagIdList,flagIdStrings,setAsActive=None):
        
        fieldname = 'Triangle Flag'
        
        # Save these for future need.
        self.flagIdList    = flagIdList
        self.flagIdStrings = flagIdStrings
        
        if setAsActive is None:
            setAsActive =  len(np.unique(flag)) > 1
            
        self.addScalarField(fieldname, flag, setAsActive=setAsActive)

    
    def addIndex(self,setAsActive=True):
        
        fieldname = 'Triangle index'
        data=np.arange(0,self.nTriangles,1,dtype=float)
        self.addScalarField(fieldname, data, setAsActive=setAsActive)
    
    def addScalarField(self,fieldname,data,setAsActive=True):
        '''
        @param fieldname: A string describing the data
        @param data: length-n array of values for coloring the triangles 
        '''

        print('Creating a VTK array representing "'+fieldname+'" and attaching it to the VTK PolyData.')
        
     
        
        field=vtk.vtkDoubleArray()
        field.SetNumberOfTuples(self.nTriangles)
        field.SetNumberOfComponents(1)
        field.SetName(fieldname)
        for i,l in enumerate(data):
            field.SetTuple(i,[l])
            
        #self.vtkPolyData.GetCellData().SetScalars(field) # This one replaces.
        self.vtkPolyData.GetCellData().AddArray(field)
        if setAsActive:
            self.vtkPolyData.GetCellData().SetActiveScalars(fieldname);
        self.vtkPolyData.Modified()
        
        
        
    def writeVtp(self, filename):
        print('Writing "' +filename+'".')
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.vtkPolyData)
        writer.Write()
        


    def readVtp(self, filename):
        print('Reading "'+filename+'".')
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        self.vtkPolyData = reader.GetOutput()


    def plot(self,manual_range=None,logarithmicColorScale=True,camControl=None,
             invertColormap=True,replaceBottom=True):
        '''
        @param manual_range: (opt) length two vector/tuple defining the coloraxis range as 10**manual_range[0] -- 10**manual_range[0] 
        @param camControl: (opt) an instance of a5VtkWall.camControl. 
               It will be associated with the new camera in the plot. The requested camera settings should be set before calling this routine.
        '''

        output = self.vtkPolyData

        # GetScalars() returns the active scalar field.
        title = output.GetCellData().GetScalars().GetName()
        print('Plotting "'+title+'".')

        scalar_range = output.GetScalarRange()

        print('Original color value range {} - {}'.format(scalar_range[0],scalar_range[1]))
        if manual_range is not None:
            #manual_range=(2,7)
            scalar_range = (10.0**manual_range[0], 10.0**manual_range[1])
            #print(scalar_range)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(output)

        mapper.SetScalarRange(scalar_range)


        '''
        Get the vtkPolyData -object. This could be useful for e.g.
        as input for
            A5vtkw = a5VtkWall(W)
            Mapper = vtk.vtkPolyDataMapper()
            Mapper.SetInputData(A5vtkw.getVtkPolyData())
        '''

        lut = vtk.vtkLookupTable()
        if logarithmicColorScale:
            lut.SetScaleToLog10()
        lut.Build()
        # Invert the colormap
        if invertColormap:
            nColors=lut.GetNumberOfTableValues()
            colors=list( range( nColors ) )
            for iColor in range( nColors ):
                colors[iColor]=lut.GetTableValue(iColor)
            for iColor in range( nColors ):
                lut.SetTableValue(nColors-iColor-1, colors[iColor] )
        # replace bottom with white
        if replaceBottom:
            lut.SetTableValue(0,(1.0,1.0,1.0,1.0))

        mapper.SetLookupTable(lut)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        ren = vtk.vtkRenderer()

        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)

        ren.GradientBackgroundOn();
        ren.SetBackground(1,1,1);
        ren.SetBackground2(0,0,0);

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        ren.AddActor(actor)

        # create the scalar_bar
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetOrientationToHorizontal()
        scalar_bar.SetLookupTable(lut)
        scalar_bar.SetTitle(title)
        if manual_range is not None:
            nTicks = int(round(manual_range[1]-manual_range[0]+1))
            if nTicks >=3:
                scalar_bar.SetNumberOfLabels(nTicks)

        # Enable user interface interactor
        print("Initializing vtkRenderWindowInteractor.")
        iren.Initialize()
        print("Rendering the vtkRenderWindow (1)")
        renWin.Render()


        # create the scalar_bar_widget
        scalar_bar_widget = vtk.vtkScalarBarWidget()
        scalar_bar_widget.SetInteractor(iren)
        scalar_bar_widget.SetScalarBarActor(scalar_bar)
        scalar_bar_widget.On()



        # Different interactor styles
        # Type "j", "t", "c" or "a" to switch
        # Or "m"??
        style = vtk.vtkInteractorStyleSwitch()
        style.SetCurrentStyleToTrackballCamera()
        iren.SetInteractorStyle(style)

        print("Rendering the vtkRenderWindow (1)")
        renWin.Render()


        ren.ResetCamera()
        if camControl is None:
            print('Setting default camera position.')
            ren.GetActiveCamera().SetPosition(6.0, 0.0, 0.0)
            ren.GetActiveCamera().SetViewUp(0.0, 0.0, 1.0)
        else:
            # Set the camera using the parameters in the camControl object given as a parameter.
            camControl.vtkCamera = ren.GetActiveCamera()
            camControl.applyValuesToVtkCamera()
            
        print("Rendering the vtkRenderWindow (2)")
        renWin.Render()


        print(self.plotHelpText)

        print("Rendering the vtkRenderWindow (3)")
        renWin.Render()

        #print("Before Start()")
        iren.Start()
        #print("After Start()")





    plotHelpText="""
    For VTK 5.0:

    vtkInteractorStyle implements the "joystick" style of interaction. That is, holding down the mouse keys generates a stream of events that cause continuous actions (e.g., rotate, translate, pan, zoom). (The class vtkInteractorStyleTrackball implements a grab and move style.) The event bindings for this class include the following:

    Keypress j / Keypress t: toggle between joystick (position sensitive) and trackball (motion sensitive) styles. In joystick style, motion occurs continuously as long as a mouse button is pressed. In trackball style, motion occurs when the mouse button is pressed and the mouse pointer moves.
    Keypress c / Keypress a: toggle between camera and actor modes. In camera mode, mouse events affect the camera position and focal point. In actor mode, mouse events affect the actor that is under the mouse pointer.
    Button 1: rotate the camera around its focal point (if camera mode) or rotate the actor around its origin (if actor mode). The rotation is in the direction defined from the center of the renderer's viewport towards the mouse position. In joystick mode, the magnitude of the rotation is determined by the distance the mouse is from the center of the render window.
    Button 2: pan the camera (if camera mode) or translate the actor (if actor mode). In joystick mode, the direction of pan or translation is from the center of the viewport towards the mouse position. In trackball mode, the direction of motion is the direction the mouse moves. (Note: with 2-button mice, pan is defined as <shift>-Button 1.)
    Button 3: zoom the camera (if camera mode) or scale the actor (if actor mode). Zoom in/increase scale if the mouse position is in the top half of the viewport; zoom out/decrease scale if the mouse position is in the bottom half. In joystick mode, the amount of zoom is controlled by the distance of the mouse pointer from the horizontal centerline of the window.
    Keypress 3: toggle the render window into and out of stereo mode. By default, red-blue stereo pairs are created. Some systems support Crystal Eyes LCD stereo glasses; you have to invoke SetStereoTypeToCrystalEyes() on the rendering window.
    Keypress e: exit the application.
    Keypress f: fly to the picked point
    Keypress p: perform a pick operation. The render window interactor has an internal instance of vtkCellPicker that it uses to pick.
    Keypress r: reset the camera view along the current view direction. Centers the actors and moves the camera so that all actors are visible.
    Keypress s: modify the representation of all actors so that they are surfaces.
    Keypress u: invoke the user-defined function. Typically, this keypress will bring up an interactor that you can type commands in.
    Keypress w: modify the representation of all actors so that they are wireframe.

    """
    # Above is from https://vtk.org/doc/release/5.0/html/a01650.html

    # For vtk 5.8, check the subpages via "see also" in https://vtk.org/doc/release/5.8/html/a01087.html



class camControl():
    
    
    def __init__(self,vtkCamera=None):
        self.focalPoint    = None  # (x,y,z)
        self.position      = None  # (x,y,z)
        self.viewUp        = None  # (x,y,z)
        self.__viewAngle     = None  # degrees (vertical)
        self.__clipRange     = None  # (min,max)
        
        self.vtkCamera = vtkCamera
        if self.vtkCamera is not None:
            self.loadValuesFromVtkCamera()

    def loadValuesFromVtkCamera(self):
        '''
        Load the settings from the VtkCamera
        '''
        if self.vtkCamera is None:
            raise ValueError('No camera to load from.')
            
        self.focalPoint = self.vtkCamera.GetFocalPoint()
        self.position   = self.vtkCamera.GetPosition()
        self.viewUp     = self.vtkCamera.GetViewUp()
        self.viewAngle  = self.vtkCamera.GetViewAngle()
    
    def applyValuesToVtkCamera(self):
        if self.vtkCamera is None:
            raise ValueError('No camera to apply to.')
        
        #self.printValues()
        self.vtkCamera.SetViewAngle( self.viewAngle )
        self.vtkCamera.SetPosition(  self.position  )
        self.vtkCamera.SetViewUp(    self.viewUp    )
        self.vtkCamera.SetFocalPoint(self.focalPoint)

        self.vtkCamera.OrthogonalizeViewUp()
        #minClip = self.vtkCamera.GetDistance() * 1.0e-3
        #maxClip = 1.0e2
        #self.vtkCamera.SetClippingRange(minClip,maxClip)
        self.vtkCamera.SetClippingRange(self.clipRange)
        #self.printValues()

    def printValues(self):
        print('Cam: '+str(self.vtkCamera))
        print('Cam: position   ({:6.3f},{:6.3f},{:6.3f})\n'.format(
            self.position[0],self.position[1],self.position[2]))
        print('Cam: focalpoint ({:6.3f},{:6.3f},{:6.3f})\n'.format(
            self.focalPoint[0],self.focalPoint[1],self.focalPoint[2]))
        print('Cam: viewUp     ({:6.3f},{:6.3f},{:6.3f})\n'.format(
            self.viewUp[0],self.viewUp[1],self.viewUp[2]))
        print('Cam: viewAngle   {:6.3f}\n'.format(
            self.viewAngle))
    
    @property
    def clipRange(self):
        return self.__clipRange

    def clipRange(self,cr):
        if len(cr) != 2:
            raise ValueError('ClipRange should be length two tuple.')
        if cr[0] >= cr[1]:
            raise ValueError('ClipRange should be increasing')
        self.__clipRange = cr

     
#     @property
#     def focalPoint(self):
#         return self.__focalPoint
#     
#     @property
#     def position(self):
#         return self.__position
#     
#     @property
#     def viewUp(self):
#         return self.__viewUp
#     
    @property
    def viewAngle(self):
        return self.__viewAngle
    
    @viewAngle.setter
    def viewAngle(self,angle):
        if angle <= 0.0:
            raise ValueError('View angle must be positive "{}".'.angle)        
        self.__viewAngle = angle
        


