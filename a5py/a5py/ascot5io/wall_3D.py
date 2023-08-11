"""
3D wall IO.

File: wall_3D.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

import a5py.wall.plot as plot

def write_hdf5(fn, nelements, x1x2x3, y1y2y3, z1z2z3, desc=None,
               flag=None, flagIdList=None, flagIdStrings=None, n=None):
    """
    Write 3D wall input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        nelements : int <br>
            Number of wall triangles
        x1x2x3 : array_like (nelements,3) <br>
            Each triangle's vertices' x coordinates [m].
        y1y2y3 : array_like (nelements,3) <br>
            Each triangle's vertices' y coordinates [m].
        z1z2z3 : array_like (nelements,3) <br>
            Each triangle's vertices' z coordinates [m].
        desc : str, optional <br>
            Input description.
        flag : array_like (nelements,1) <br>
            Integer array depicting the wall component of each triangle,
        flagIdList, array_like (nUniqueFlags) <br>
            List of keys of the flagIdStrings (int)
        flagIdStrings, array_like (nUniqueFlags) <br>
            List of values of the flagIdStrings (string)
        n, optional, ignored

    Returns:
        Name of the new input that was written.
    """

    assert x1x2x3.shape == (nelements,3)
    assert y1y2y3.shape == (nelements,3)
    assert z1z2z3.shape == (nelements,3)

    # Replace what ever we happen to get with a plain integer.
    # Seems to be some array sometimes
    nelements = x1x2x3.shape[0]

    if flag is None:
        flag = np.zeros(shape=(nelements,1),dtype=int)
    else:
        assert ( flag.shape == (nelements,1) or flag.shape == (nelements,) )

    parent = "wall"
    group  = "wall_3D"
    gname  = ""

    # Convert strings to the favorite format
    if flagIdStrings is not None:
        strlen = 0
        for s in flagIdStrings:
            if len(s) >  strlen:
                strlen = len(s)
        fids=np.empty(shape=(len(flagIdStrings),),dtype='|S{}'.format(strlen) )
        for istr,s in enumerate(flagIdStrings):
            fids[istr]=s

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset('x1x2x3',     (nelements,3), data=x1x2x3,    dtype='f8')
        g.create_dataset('y1y2y3',     (nelements,3), data=y1y2y3,    dtype='f8')
        g.create_dataset('z1z2z3',     (nelements,3), data=z1z2z3,    dtype='f8')
        g.create_dataset('nelements',  (1,1),         data=nelements, dtype='i4')
        fl = g.create_dataset('flag',(nelements,1),   data=np.reshape(flag,(nelements,1)),      dtype='i4')

        if flagIdList is not None and flagIdStrings is not None:
            fl.attrs.create(name='flagIdList',    data=np.array( flagIdList    ),  dtype='i4')
            fl.attrs.create(name='flagIdStrings', data=fids, dtype=h5py.string_dtype(encoding='utf-8',length=None)           )



    return gname


def read_hdf5(fn, qid):
    """
    Read 3D wall input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "wall/wall_3D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

        nTriangles = out['x1x2x3'].shape[0]

        if path+'/flag' in f:
            flagAttrs = ['flagIdStrings','flagIdList']
            for s in flagAttrs:
                if s in f[ path+'/flag' ].attrs:
                    out[s] = f[ path+'/flag' ].attrs.get(s)

        if not 'n' in out:
            out['n'] = np.array([nTriangles])

        if not 'flag' in out:
            out['flag'] = np.zeros(shape=(nTriangles,),dtype=int)

        if 'flagIdStrings' in out:
            # We need to decode the bytearrays into strings.
            s=[]
            for S in out['flagIdStrings']:
                if not isinstance(S,str):
                    s.append(S.decode('utf-8'))
                else:
                    s.append(S)
            out['flagIdStrings']=s
        else:
            # Generate some flag names.
            out['flagIdList']=np.unique(out['flag'])
            out['flagIdStrings']=[]
            for fl in out['flagIdList']:
                out['flagIdStrings'].append('Flag {}'.format(fl))

    return out


class wall_3D(AscotData):
    """
    Object representing wall_3D data.
    """

    number_of_elements = None

    def read(self):
        W=read_hdf5(self._file, self.get_qid())

        # The dimensionality of  nTriangles = a5wall['nelements'][0][0] varies, so by-pass it.
        self.number_of_elements = W['x1x2x3'].shape[0]
        return W

    def getNumberOfElements(self):
        if self.number_of_elements is None:
            self.read()
        return self.number_of_elements


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    def remove_small_triangles(self,maximumAreaToRemove=0.0,data=None):
        "The modification happens in-place. No deep copy is made!"

        if data is None:
            w = self.read()
        else:
            w = data

        A = self.area(data=w)

        keep =  ( A > maximumAreaToRemove )

        fields = ["x1x2x3","y1y2y3","z1z2z3"]
        for f in fields:
            w[f] = w[f][keep,:]
        w['flag'] = w['flag'][keep]

        nOld = w['n'][:]
        nNew = len(w['flag'])
        w['n'][:]         = nNew
        w['nelements'][:] = nNew

        print('Removing {}/{} triangles'.format(nOld-nNew,nOld))

        return w

    def area(self, data= None):
        if data is None:
            w = self.read()
        else:
            w = data


        ab_x = w["x1x2x3"][:,1] - w["x1x2x3"][:,0]
        ab_y = w["y1y2y3"][:,1] - w["y1y2y3"][:,0]
        ab_z = w["z1z2z3"][:,1] - w["z1z2z3"][:,0]

        ac_x = w["x1x2x3"][:,2] - w["x1x2x3"][:,0]
        ac_y = w["y1y2y3"][:,2] - w["y1y2y3"][:,0]
        ac_z = w["z1z2z3"][:,2] - w["z1z2z3"][:,0]

        A = 0.5 * np.sqrt(   (ab_y * ac_z - ab_z * ac_y)**2
                           + (ab_z * ac_x - ab_x * ac_z)**2
                           + (ab_x * ac_y - ab_y * ac_x)**2 )

        return A

    def getAspointsAndVertices(self, removeDuplcatePoints=True):
        a5wall = self.read()
        nTriangles = self.number_of_elements

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

        if not removeDuplcatePoints:
            return points,vertices

        points,vertices = self._removeDuplicatePoints(points,vertices)
        return points,vertices

    def _removeDuplicatePoints(self,points,vertices):
        ''' Remove duplicate points from points and vertices representation of the 3D wall
        '''

        upoints,inverse = np.unique(points, return_index=False,return_inverse=True, axis=0)


        uvertices = inverse[vertices]

        return upoints,uvertices


    def plotRz(self, axes=None, phi=None):
        data = self.read()
        if phi is not None:
            plot.plot_intersection(data["x1x2x3"], data["y1y2y3"],
                                   data["z1z2z3"], phi, axes=axes)
        else:
            plot.plot_projection(data["x1x2x3"], data["y1y2y3"],
                                 data["z1z2z3"], axes=axes)


    def toVtk(self):
        import a5py.wall.a5vtkwall

        a5VTKwall = a5py.wall.a5vtkwall.a5VtkWall()

        a5VTKwall.fromA5wall(self)

        W=self.read()

        # Add index of running triangle number
        a5VTKwall.addIndex(setAsActive=True)

        # Add coloring based on flags. (Get's activated only if there are more than one flag.
        a5VTKwall.addFlag(W['flag'],W['flagIdList'],W['flagIdStrings'])

        return a5VTKwall

    def append(self,newWall,data=None):
        '''
        Returns the wall  with newWall["x1x2x3], newWall["y1y2y3"],  newWall["z1z2z3"], newWall["flag"] to the existing triangles.
        If data=None (or not given) reads the old triangles from file.
        '''
        if data is None:
            w = self.read()
        else:
            w = data

        fields = ["x1x2x3","y1y2y3","z1z2z3","flag"]
        for f in fields:
            #print('Field "{}" has the shapes (data,newWall):'.format(f))
            #print(w[f].shape, newWall[f].shape)


            if(len(w[f].shape)==1):
                if( len(newWall[f]) != len(newWall[f].ravel()) ):
                    raise Exception("Wrong size field {}".format(f))
                w[f] = np.concatenate( (w[f], newWall[f].ravel() ), axis=0 )
            else:
                w[f] = np.concatenate( (w[f], newWall[f]         ), axis=0 )



        w['n'] = w['n'] + newWall['n']

        return w


    def move_component(self, movement, direction, component):
        '''
        @Params:
        filename    Filename where the wall will be saved.
        movement    How much the component is moved in metres
        direction    Which direction the component is moved [x, y, z]
        component     Bool array of the points of a component
        wall         Wall object where the wall is copied from, if left empty active
        wall from filename will be used

        @return the new wall data
        '''


        rwall = self.read()


        #Movement in metres
        t = movement/np.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2)

        rwall['x1x2x3'][component,:] += t*direction[0]
        rwall['y1y2y3'][component,:] += t*direction[1]
        rwall['z1z2z3'][component,:] += t*direction[2]

        new_rwall = {'x1x2x3':rwall['x1x2x3'],'y1y2y3':rwall['y1y2y3'],'z1z2z3':rwall['z1z2z3'],
                'desc':'ICRHmv{}m'.format(movement),
            'flag':rwall['flag'], 'flagIdList':rwall['flagIdList'],
            'flagIdStrings':rwall['flagIdStrings'], 'nelements':rwall['n']}

        return new_rwall


    def set_flags( self, names, indexes):
        '''
        @Params:
        names       list of strings containing names of various flags
        indexes     which flag indexes (zero indexed, last index value not included) correspond to each flag, list of tuples

        example with 1000 triangles:
                 names   = ['foo','bar','baz']
                 indexes = [(0,10), (10,100), (100,1000)]

        @return the new wall data
        '''
        W = self.read()

        flagIdList    = []
        flagIdStrings = []

        for i,name in enumerate(names):
            W['flag'][indexes[i][0]:(indexes[i][1])] = i
            flagIdList.append(i)
            flagIdStrings.append(name)

        W['flagIdList']    = np.array(  flagIdList )
        W['flagIdStrings'] = flagIdStrings

        return W

