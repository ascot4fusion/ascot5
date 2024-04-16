"""Input representing first wall and other plasma facing components.

Wall input is used in simulations to check intersections between a marker
trajectory and physical components like first wall of fast-ion-loss detector.
It is also used in BBNBI5 to estimate shinethrough and in ASCOT-BMC.

Wall models can also be used to mimic LCFS to catch markers that escape plasma.
This could prove more detailed information than just using the MAXRHO end
condition.
"""
import h5py
import numpy as np
import unyt

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

import a5py.physlib as physlib
from a5py.routines.plotting import pv

class wall_2D(DataGroup):
    """Contour in Rz-plane that represents an axisymmetric wall.

    This wall is simple and fast in simulations, but it does not allow
    evaluation of wall loads. For those, use 3D wall instead.

    The wall doesn't have to form a closed loop but it prevents markers
    from escaping the computational domain.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]
                if key == "nelements":
                    out[key] = int(out[key])

        return out

    @staticmethod
    def write_hdf5(fn, nelements, r, z, desc=None):
        """Write input data to the HDF5 file.

        First vertex shouldn't correspond to the last vertex, i.e., don't give
        (R,z) coordinates that make a closed loop.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        nelements : int
            Number of wall segments.
        r : array_like (nelements,1)
            R coordinates of wall segment vertices [m].
        z : array_like (nelements,1)
            z coordinates of wall segment vertices [m].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If the input is inconsistent.
        """
        if r.size != nelements:
            raise ValueError("Size of r does not match the number of elements.")
        if z.size != nelements:
            raise ValueError("Size of z does not match the number of elements.")

        parent = "wall"
        group  = "wall_2D"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("r",         (nelements,1), data=r, dtype='f8')
            g.create_dataset("z",         (nelements,1), data=z, dtype='f8')
            g.create_dataset("nelements", (1,1),         data=nelements,
                             dtype='i4')

        return gname

    def getwallcontour(self, phi=0*unyt.deg):
        """Return a cross section of the wall with a given poloidal plane.

        Parameters
        ----------
        phi : float
            Toroidal angle of the plane.

        Returns
        -------
        lines : array_like (n,2,2)
            Line segments [[[r1,z1], [r2,z2]], ...] that form the cross section.
        """
        w = self.read()
        r = np.append(w["r"], w["r"][0])
        z = np.append(w["z"], w["z"][0])
        lines = np.zeros((r.size-1, 2, 2))
        for i in range(r.size-1):
            lines[i, 0, 0] = r[i]
            lines[i, 0, 1] = z[i]
            lines[i, 1, 0] = r[i+1]
            lines[i, 1, 1] = z[i+1]
        return lines

    def getwalloutline(self, phigrid):
        """Return minimum and maximum wall R-coordinates as a function of phi.

        Parameters
        ----------
        phigrid : array_like
            Toroidal coordinates.

        Returns
        -------
        rmin : array_like, (nphi,)
            Minimum R-coordinate at each toroidal angle.
        rmax : array_like, (nphi,)
            Maximum R-coordinate at each toroidal angle.
        """
        w = self.read()
        rmin = np.amax(w["r"]) * np.ones(phigrid.shape).ravel()
        rmax = np.amin(w["r"]) * np.ones(phigrid.shape).ravel()

        return rmin, rmax

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is a very large rectangular wall.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        r = np.array([0.01, 100, 100, 0.01])
        z = np.array([-100, -100, 100, 100])
        return {"nelements":4, "r":r, "z":z,}


class wall_3D(DataGroup):
    """Triangle mesh that represents a 3D wall.

    This wall model can be obtained easily from a STEP file or from a CAD model
    with a little effort. If 3D model of the wall is not available, consider
    using the mock-up 3D wall generated by "rotating" the 2D wall contour
    toroidally in steps.

    3D wall model allows estimating the wall loads, but keep in mind that
    smaller triangles means more markers are required in simulation to gain
    sufficient statistics. On the other hand, too coarse model means that
    important details on wall load distribution could be missed.

    This wall doesn't have to form a closed "waterproof" structure, and one
    can use this to model FILD.
    """

    number_of_elements = None

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]
                if key == "nelements":
                    out[key] = int(out[key])

            nTriangles = out['x1x2x3'].shape[0]
            self.number_of_elements = nTriangles

            if path+'/flag' in f:
                flagAttrs = ['flagIdStrings','flagIdList']
                for s in flagAttrs:
                    if s in f[ path+'/flag' ].attrs:
                        out[s] = f[ path+'/flag' ].attrs.get(s)

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

    def getNumberOfElements(self):
        if self.number_of_elements is None:
            self.read()
        return self.number_of_elements

    def getwalloutline(self, z=0):
        """Return minimum and maximum wall R-coordinates as a function of phi
        at the given toroidal plane.

        Parameters
        ----------
        phigrid : array_like
            Toroidal coordinates.

        Returns
        -------
        rmin : array_like, (nphi,)
            Minimum R-coordinate at each toroidal angle.
        rmax : array_like, (nphi,)
            Maximum R-coordinate at each toroidal angle.
        """
        s1 = self.tomesh()
        xmin = np.amin(s1.points[:,0])
        xmax = np.amax(s1.points[:,0])
        ymin = np.amin(s1.points[:,1])
        ymax = np.amax(s1.points[:,1])
        r0 = 0.0
        z0 = 0.0
        dx = xmax - xmin
        dy = ymax - ymin

        s2 = pv.Plane(center=(0,0,z), direction=(0,0,1),
                      i_size=dx*1.1, j_size=dy*1.1,
                      i_resolution=1, j_resolution=1).triangulate()

        cut,_,_ = s2.intersection(s1, split_first=False, split_second=False)
        i0 = cut.lines[1::3]
        i1 = cut.lines[2::3]
        lines = np.array( [
            ( (cut.points[i0[i],0], cut.points[i0[i],1]),
              (cut.points[i1[i],0], cut.points[i1[i],1]) ) \
            for i in range(i0.size) ] )
        return lines

    def noderepresentation(self, w_indices=None):
        """
        Return an array of vertices and indices that define this wall.

        Parameters
        ----------
        w_indices : array_like, optional
            list of of triangle indices for which the node
            representation is done i.e. [ind1, ind2, ..., indn].

            Can also be a Boolean array with the same size as the number of
            wall triangles.

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
            if w_indices is None:
                w_indices = np.s_[:]

            x1x2x3 = h5["x1x2x3"][w_indices, :]
            y1y2y3 = h5["y1y2y3"][w_indices, :]
            z1z2z3 = h5["z1z2z3"][w_indices, :]
            ntriangle = x1x2x3.shape[0]

            faces     = np.zeros((ntriangle,4), dtype="i8")
            vertices  = np.zeros((ntriangle*3,3), dtype="f8")

            faces[:, 0] = 3
            faces[:, 1] = np.arange(ntriangle) * 3 + 0
            faces[:, 2] = np.arange(ntriangle) * 3 + 1
            faces[:, 3] = np.arange(ntriangle) * 3 + 2

            vertices[0::3,0] = x1x2x3[:,0]
            vertices[0::3,1] = y1y2y3[:,0]
            vertices[0::3,2] = z1z2z3[:,0]
            vertices[1::3,0] = x1x2x3[:,1]
            vertices[1::3,1] = y1y2y3[:,1]
            vertices[1::3,2] = z1z2z3[:,1]
            vertices[2::3,0] = x1x2x3[:,2]
            vertices[2::3,1] = y1y2y3[:,2]
            vertices[2::3,2] = z1z2z3[:,2]

        return (vertices, faces)

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

    def area(self, normal=False, data=None):
        """Calculate wall element area.

        Parameters
        ----------
        normal : bool, optional
            If True, calculate and return the surface normal vectors as well.
        data : dict, optional
            Dictionary with the wall data. If ``None``, the data is read from
            the file.

        Returns
        -------
        area : array_like, (nelement,)
            Surface areas of the wall elements.
        nvec : array_like, (nelement,3)
            Normal vectors of the wall elements in cartesian basis.

            Note that the direction of the normal vector is not specifiec (it
            can vary between elements and point either inside or outside).
        """
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

        area = 0.5 * np.sqrt(   (ab_y * ac_z - ab_z * ac_y)**2
                              + (ab_z * ac_x - ab_x * ac_z)**2
                              + (ab_x * ac_y - ab_y * ac_x)**2 ) * unyt.m**2
        if not normal: return area

        nvec = np.array([
            ab_y * ac_z - ab_z * ac_y,
            ab_z * ac_x - ab_x * ac_z,
            ab_x * ac_y - ab_y * ac_x
        ])
        nvec /= np.sqrt(np.sum(nvec**2, axis=0))
        return area, nvec

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

    def getwallcontour(self, phi=0*unyt.deg):
        """Return a cross section of the wall with a given poloidal plane.

        Parameters
        ----------
        phi : float
            Toroidal angle of the plane.

        Returns
        -------
        lines : array_like (n,2,2)
            Line segments [[[r1,z1], [r2,z2]], ...] that form the cross section.
        """
        phi = phi.to("rad").v
        s1 = self.tomesh()
        rmin = np.amin(np.sqrt( s1.points[:,0]**2 + s1.points[:,1]**2 ))
        rmax = np.amax(np.sqrt( s1.points[:,0]**2 + s1.points[:,1]**2 ))
        zmin = np.amin(s1.points[:,2])
        zmax = np.amax(s1.points[:,2])
        r0 = 0.5*(rmin+rmax)
        z0 = 0.5*(zmin+zmax)
        dr = (rmax - rmin)
        dz = (zmax - zmin)

        x, y, z = physlib.pol2cart(r0, phi, z0)
        s2 = pv.Plane(center=(x,y,z), direction=(-np.sin(phi),np.cos(phi),0),
                      i_size=dz*1.1, j_size=dr*1.1,
                      i_resolution=1, j_resolution=1).triangulate()

        cut,_,_ = s2.intersection(s1, split_first=False, split_second=False)
        r, phi, z = physlib.cart2pol(cut.points[:,0], cut.points[:,1],
                                     cut.points[:,2])
        i0 = cut.lines[1::3]
        i1 = cut.lines[2::3]
        lines = np.array( [
            ( (r[i0[i]], z[i0[i]]), (r[i1[i]], z[i1[i]]) ) \
            for i in range(i0.size) ] )
        return lines

    def tomesh(self):
        import pyvista as pv
        return pv.PolyData( *self.noderepresentation() )

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


    @staticmethod
    def write_hdf5(fn, nelements, x1x2x3, y1y2y3, z1z2z3, flag=None,
                   flagIdList=None, flagIdStrings=None, desc=None):
        """Write input data to the HDF5 file.

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
        flag : array_like (nelements,1), optional
            Integer array depicting the wall component of each triangle.
        flagIdList : array_like (nUniqueFlags), optional
            List of keys (int) of the flagIdStrings.
        flagIdStrings : array_like (nUniqueFlags), optional
            List of values (str) of the flagIdStrings.
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If the triangle vertices or flags have incorrect shape.
        """
        if type(nelements) == int:
            pass
        elif type(nelements) == list:
            while type(nelements) == list:
                nelements = nelements[0]
        elif type(nelements) == np.ndarray:
            nelements = int(nelements)
        else:
            raise ValueError("Unsupported format for nelements.")

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
            flag = np.zeros(shape=(nelements,),dtype=int)
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
            fids = np.empty( shape=(len(flagIdStrings),),
                             dtype='|S{}'.format(strlen) )
            for istr,s in enumerate(flagIdStrings):
                fids[istr]=s

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nelements', (1,1), data=nelements, dtype='i4')
            g.create_dataset('x1x2x3', (nelements,3), data=x1x2x3, dtype='f8')
            g.create_dataset('y1y2y3', (nelements,3), data=y1y2y3, dtype='f8')
            g.create_dataset('z1z2z3', (nelements,3), data=z1z2z3, dtype='f8')

            fl = g.create_dataset('flag', (nelements,1), data=flag, dtype='i4')
            if flagIdList is not None and flagIdStrings is not None:
                flagIdList = np.array(flagIdList)
                fl.attrs.create(name='flagIdList', data=flagIdList, dtype='i4')
                fl.attrs.create(name='flagIdStrings', data=fids,
                                dtype=h5py.string_dtype(encoding='utf-8',
                                                        length=None))

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output consists of a single triangle on a poloidal plane.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return {"nelements":1, "x1x2x3":np.array([[4],[8],[8]]).T,
                "y1y2y3":np.array([[0],[0],[0]]).T,
                "z1z2z3":np.array([[-4],[-4],[4]]).T}

    @staticmethod
    def convert_wall_2D(nphi, **kwargs):
        """Convert :class:`wall_2D` input to :class:`wall_3D` input.

        This function repeats the 2D wall contour around the torus at intervals
        Delta phi = 2 pi / (nphi - 1). Then each (R, z) point is connected to a
        equivalent point in the adjacent contours by a line segment forming 3D
        wall made of rectangles. Finally each rectangle is divided into two
        triangles.

        Parameters
        ----------
        nphi : int
            Number of toroidal segments to be created.
        **kwargs
            Arguments passed to :meth:`wall_2D.write_hdf5` excluding ``fn`` and
            ``desc``.

        Returns
        -------
        out : dict
            :class:`wall_2D` converted as an input for :meth:`write_hdf5`.
        """
        r = kwargs["r"]
        z = kwargs["z"]
        n = kwargs["nelements"]
        x1x2x3 = np.ones((2*n*nphi, 3))
        y1y2y3 = np.ones((2*n*nphi, 3))
        r1r2r3 = np.ones((2*n*nphi, 3))
        p1p2p3 = np.ones((2*n*nphi, 3))
        z1z2z3 = np.ones((2*n*nphi, 3))

        def pol2cart(rho, phi):
            """Poloidal coordinates to Cartesian.
            """
            x = rho * np.cos(phi)
            y = rho * np.sin(phi)
            return(x, y)

        pv = np.linspace(0, 2*np.pi, nphi+1)
        for i in range(1,nphi+1):
            for j in range(n):
                r1r2r3[(i-1)*2*n + 2*(j-1),:] = [ r[j],    r[j],  r[j-1] ]
                p1p2p3[(i-1)*2*n + 2*(j-1),:] = [ pv[i-1], pv[i], pv[i] ]
                z1z2z3[(i-1)*2*n + 2*(j-1),:] = [ z[j],    z[j],  z[j-1] ]

                r1r2r3[(i-1)*2*n + 2*(j-1) + 1,:] = [ r[j],    r[j-1],  r[j-1] ]
                p1p2p3[(i-1)*2*n + 2*(j-1) + 1,:] = [ pv[i-1], pv[i-1], pv[i] ]
                z1z2z3[(i-1)*2*n + 2*(j-1) + 1,:] = [ z[j],    z[j-1],  z[j-1] ]

        x1x2x3, y1y2y3 = pol2cart(r1r2r3, p1p2p3)

        return {"nelements" : 2*n*nphi, "x1x2x3" : x1x2x3,
                "y1y2y3" : y1y2y3, "z1z2z3" : z1z2z3}
