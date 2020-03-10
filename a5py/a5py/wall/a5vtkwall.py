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
        @param W: (optional) The output of  a5py.ascot5io.wall_3D.read_hdf5()
        No reduction of points.
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

    def addIndex(self):
        
        fieldname = 'Triangle index'
        data=np.arange(0,self.nTriangles,1,dtype=float)
        self.addScalarField(fieldname, data)
    
    def addScalarField(self,fieldname,data):
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
