'''
Created on Mar 10, 2020

@author: sjjamsa
'''

def save3DAsVTP(hdf5,filename):
    '''
    Calculate and save the wall 3D loads to a VTK file. 
    @param hdf5: a5py.ascot5io.ascot5.Ascot(filename)
    @param filename: The output filename, usually with the suffix .vtp
    '''
    
    import a5py.wallloads.calculate
    import a5py.wall.a5vtkwall
    
    wallLoad = a5py.wallloads.calculate.wallLoad3DEndstate(hdf5)
    wall = hdf5.active.wall
    
    avtk = a5py.wall.a5vtkwall.a5VtkWall(wall) 
    avtk.addIndex() # This is just the triangle index, included for just in case...
    
    fieldname = 'Wall load (W/m^2)'
    avtk.addScalarField(fieldname=fieldname,data=wallLoad)
    
    avtk.writeVtp(filename)