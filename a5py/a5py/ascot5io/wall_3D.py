"""3D wall IO.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

import a5py.wall.plot as plot


class wall_3D(DataGroup):
    """
    Object representing wall_3D data.
    """

    number_of_elements = None

    def read(self):
        W=read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

        # The dimensionality of  nTriangles = a5wall['nelements'][0][0] varies, so by-pass it.
        self.number_of_elements = W['x1x2x3'].shape[0]
        return W

    def getNumberOfElements(self):
        if self.number_of_elements is None:
            self.read()
        return self.number_of_elements


    def noderepresentation(self):
        """
        Return an array of vertices and indices that define this wall.

        Returns:
            vertices : float, array_like <br>
                Triangle vertices in a format [[x1, y1, z1], ..., [xn, yn, zn]]
                where n is 3 x number of triangles.
            faces : int, array_like <br>
                Array of dimension (number of triangles, 4) where each four
                length array defines a wall elemenet as [number of vertices, i1,
                i2, i3] where number of vertices is always 3 and i1, i2, and i3
                define the indices of the vertices in the vertices array.
        """
        with self as h5:
            ntriangle = int(h5["nelements"][:])
            faces     = np.zeros((ntriangle,4), dtype="i8")
            vertices  = np.zeros((ntriangle*3,3), dtype="f8")

            faces[:, 0] = 3
            faces[:, 1] = np.arange(ntriangle) * 3 + 0
            faces[:, 2] = np.arange(ntriangle) * 3 + 1
            faces[:, 3] = np.arange(ntriangle) * 3 + 2

            vertices[0::3,0] = h5["x1x2x3"][:,0]
            vertices[0::3,1] = h5["y1y2y3"][:,0]
            vertices[0::3,2] = h5["z1z2z3"][:,0]
            vertices[1::3,0] = h5["x1x2x3"][:,1]
            vertices[1::3,1] = h5["y1y2y3"][:,1]
            vertices[1::3,2] = h5["z1z2z3"][:,1]
            vertices[2::3,0] = h5["x1x2x3"][:,2]
            vertices[2::3,1] = h5["y1y2y3"][:,2]
            vertices[2::3,2] = h5["z1z2z3"][:,2]

        return (vertices, faces)


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

    def area(self, data=None):
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
        """
        Remove duplicate points from points and vertices representation.
        """

        upoints,inverse = np.unique(points, return_index=False,
                                    return_inverse=True, axis=0)


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

        # Add coloring based on flags. (Get's activated only if there are more
        # than one flag.
        a5VTKwall.addFlag(W['flag'],W['flagIdList'],W['flagIdStrings'])

        return a5VTKwall


    def append(self,newWall,data=None):
        """
        Returns the wall with newWall["x1x2x3], newWall["y1y2y3"], newWall["z1z2z3"], newWall["flag"] to the existing triangles.
        If data=None (or not given) reads the old triangles from file.
        """
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
        """
        @Params:
        filename    Filename where the wall will be saved.
        movement    How much the component is moved in metres
        direction    Which direction the component is moved [x, y, z]
        component     Bool array of the points of a component
        wall         Wall object where the wall is copied from, if left empty active 
        wall from filename will be used

        @return the new wall data
        """
        rwall = self.read()

        #Movement in metres
        t = movement/np.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2) 

        rwall['x1x2x3'][component,:] += t*direction[0]
        rwall['y1y2y3'][component,:] += t*direction[1]
        rwall['z1z2z3'][component,:] += t*direction[2]

        new_rwall = {'x1x2x3':rwall['x1x2x3'],'y1y2y3':rwall['y1y2y3'],
                     'z1z2z3':rwall['z1z2z3'],
                     'desc':'ICRHmv{}m'.format(movement),
            'flag':rwall['flag'], 'flagIdList':rwall['flagIdList'],
            'flagIdStrings':rwall['flagIdStrings'], 'nelements':rwall['n']}

        return new_rwall

    @staticmethod
    def write_hdf5(fn, nelements, x1x2x3, y1y2y3, z1z2z3, desc=None,
               flag=None, flagIdList=None, flagIdStrings=None):
        """Write 3D wall input in HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        nelements : int
            Number of wall triangles
        x1x2x3 : array_like (nelements,3)
            Each triangle's vertices' x coordinates [m].
        y1y2y3 : array_like (nelements,3)
            Each triangle's vertices' y coordinates [m].
        z1z2z3 : array_like (nelements,3)
            Each triangle's vertices' z coordinates [m].
        desc : str, optional
            Input description.
        flag : array_like (nelements,1), optional
            Integer array depicting the wall component of each triangle.
        flagIdList : array_like (nUniqueFlags), optional
            List of keys (int) of the flagIdStrings.
        flagIdStrings : array_like (nUniqueFlags), optional
            List of values (str) of the flagIdStrings.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If the triangle vertices or flags have incorrect shape.
        """
        if x1x2x3.shape != (nelements,3):
            raise ValueError(
                "Shape of x1x2x3 was " + str(x1x2x3.shape) + " but expected ("
                + str(nelements) + ",3)")
        if y1y2y3.shape != (nelements,3):
            raise ValueError(
                "Shape of y1y2y3 was " + str(y1y2y3.shape) + " but expected ("
                + str(nelements) + ",3)")
        if z1z2z3.shape != (nelements,3):
            raise ValueError(
                "Shape of z1z2z3 was " + str(z1z2z3.shape) + " but expected ("
                + str(nelements) + ",3)")

        if flag is None:
            flag = np.zeros(shape=(nelements,1),dtype=int)
        elif flag.shape != (nelements,1):
            raise ValueError(
                "Shape of flag was " + str(flag.shape) + " but expected ("
                + str(nelements) + ",1)")

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

            g.create_dataset('x1x2x3',    (nelements,3), data=x1x2x3,    dtype='f8')
            g.create_dataset('y1y2y3',    (nelements,3), data=y1y2y3,    dtype='f8')
            g.create_dataset('z1z2z3',    (nelements,3), data=z1z2z3,    dtype='f8')
            g.create_dataset('nelements', (1,1),         data=nelements, dtype='i4')

            fl = g.create_dataset('flag', (nelements,1), data=flag,      dtype='i4')
            if flagIdList is not None and flagIdStrings is not None:
                flagIdList = np.array(flagIdList)
                fl.attrs.create(name='flagIdList', data=flagIdList, dtype='i4')
                fl.attrs.create(name='flagIdStrings', data=fids)

        return gname

    @staticmethod
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
                    s.append(S.decode('utf-8'))
                out['flagIdStrings']=s
            else:
                # Generate some flag names.
                out['flagIdList']=np.unique(out['flag'])
                out['flagIdStrings']=[]
                for fl in out['flagIdList']:
                    out['flagIdStrings'].append('Flag {}'.format(fl))

        return out
