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


    def __init__(self, a5wall=None):
        '''
        @param a5wall: (optional) The output of  a5py.ascot5io.wall_3D.read_hdf5()
        '''
        
        if a5wall is not None:
            self.fromA5wall(a5wall)
         
    
        
    def fromA5wall(self, a5wall):
        '''
        @param a5wall: (optional) The output of  a5py.ascot5io.wall_3D.read_hdf5()
        No reduction of points.
        '''
        
        points,vertices  = self.pointsAndVerticesFromA5wall(a5wall)
        
        self.fromPointsAndVertices(points,vertices)
    
    def pointsAndVerticesFromA5wall(self,a5wall):
    
        nTriangles = a5wall['nelements'][0][0]

        # Create an array 
        points = np.zeros(shape=(3*nTriangles,3))
        
        points[ ::3,0]=a5wall['x1x2x3'][:,0]
        points[1::3,0]=a5wall['x1x2x3'][:,1]
        points[2::3,0]=a5wall['x1x2x3'][:,2]
        points[ ::3,1]=a5wall['y1y2y3'][:,0]
        points[1::3,1]=a5wall['y1y2y3'][:,1]
        points[2::3,1]=a5wall['y1y2y3'][:,2]
        points[ ::3,2]=a5wall['z1z2z3'][:,0]
        points[1::3,2]=a5wall['z1z2z3'][:,1]
        points[2::3,2]=a5wall['z1z2z3'][:,2]                     

    
        vertices = np.reshape( np.arange(0,3*nTriangles,1,dtype=int), (nTriangles,3) )
        
       
        
        return points,vertices
    
    
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
