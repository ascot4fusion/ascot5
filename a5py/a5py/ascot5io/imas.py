import numpy as np

class a5imas:


    def __init__(self):
        self.ids_name = "ids" #This should be overwritten by dervied classes

    def open(self, user, tokamak, version, shot, run, occurrence=0):

        # Put this inside the function, not to disturb usage where imas is not available
        import imas

        # ./idsdump.py g2jvarje test 3 92436 0272 distributions
        # user = "g2jvarje"
        # tokamak = "test"
        # version = "3"
        # shot = 92436
        # run = 272
        # ids = "wall"



        self.ids = imas.ids(shot, run)
        self.ids.open_env(user, tokamak, version)

        #ids = ids.split('/')
        #if len(ids) == 1:
        #    occurrence = 0
        #else:
        #    occurrence = int(ids[1])
        #ids = ids[0]

        ids = self.ids_name

        idsdata = self.ids.__dict__[self.ids_name]
        if 'get' not in dir(idsdata):
            idsdata = self.ids.__dict__[self.ids_name + 'Array']
            idsdata.get(occurrence)
            #if idsdata.array:
            #    for slice in idsdata.array:
            #        print(slice)
        else:
            idsdata.get(occurrence)
            #print(idsdata)

        return self.ids

    def close(self):
        self.ids.close()


class wall_2d(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "wall"


    def read(self, user, tokamak, version, shot, run, occurrence=0 ):
        """
        Read an IMAS 2D-wall into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from
              wall.description_2d[0].limiter.unit[0].outline.r
        and
              wall.description_2d[0].limiter.unit[0].outline.z

        """

        # Put this inside the function, not to disturb usage where imas is not available
        import imas


        timeIndex = 0
        unit      = 0

        itm = self.open( user, tokamak, version, shot, run, occurrence )

        #itm.wall.get()
        r = self.ids.wall.description_2d[timeIndex].limiter.unit[unit].outline.r
        z = self.ids.wall.description_2d[timeIndex].limiter.unit[unit].outline.z

        self.close()

        w = {
            "r"         : r,
            "z"         : z,
            "nelements" : len(r),
        }

        return w

class wall_3d(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "wall"


    def read(self, user, tokamak, version, shot, run, occurrence=0 ):
        """
        Read an IMAS 3D-wall into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from
              wall.description_ggd[0].grid_ggd[0].grid_ggd.space[0]
        Flag is set to the index of description_ggd == 0

        """

        itm = self.open( user, tokamak, version, shot, run, occurrence )

        description_ggd_index = 0
        grid_ggd_index        = 0
        space_index           = 0
        
        
        desc_ggd = self.ids.wall.description_ggd[description_ggd_index]
        grid_ggd = desc_ggd.grid_ggd[grid_ggd_index]
        space    = grid_ggd.space[space_index]
        
        if not np.all( space.coordinates_type == np.array([1, 2, 3]) ):
            raise ValueError("Space coordinates not [1,2,3] in IMAS (3d)wall ids u:{} db:{} v:{} s:{} r:{} o:{} desc:{} grid:{} sp:{}".format(
                user,tokamak,version,shot,run,occurrence,description_ggd_index,grid_ggd_index,space_index ))


        node_coordinates,edge_indexes,tria_indexes = self.get_nodes_edges_polygons(space)

        self.close()


        # Convert 3D wall data to ASCOT5 format
                             
        xyz = self.get_xyz(node_coordinates,tria_indexes)
        

        nelements = xyz.shape[0]
        flag      = np.ones( (nelements,1) ) * description_ggd_index
        w = {
            "x1x2x3"         : xyz[:,:,0],
            "y1y2y3"         : xyz[:,:,1],
            "z1z2z3"         : xyz[:,:,2],
            "nelements"      : nelements,
            "flag"           : flag,
            "flagIdList"     : flag,
            "flagIdStrings"  : ["IMAS 3D-wall ids u:{} db:{} v:{} s:{} r:{} o:{} desc:{} grid:{} sp:{}".format(
                user,tokamak,version,shot,run,occurrence,description_ggd_index,grid_ggd_index,space_index )],

        }

        return w


    def get_nodes_edges_polygons(self,space):
        '''
        Copy the node coordinates, edge indexes and polygons from the space to separae numpy arrays:

        returns node_coordinates[3,n_nodes], edge_indexes[2,n_edges],tria_indexes[3,n_tris]

        '''


        n_nodes=len(space.objects_per_dimension[0].object)
        n_edges=len(space.objects_per_dimension[1].object)
        n_tris =len(space.objects_per_dimension[2].object)
        node_coordinates = np.zeros( (3,n_nodes),dtype=float )
        edge_indexes     = np.zeros( (2,n_edges),dtype=int   )
        tria_indexes     = np.zeros( (3,n_tris), dtype=int   )

        for i in range(n_nodes):
            node_coordinates[:,i] = space.objects_per_dimension[0].object[i].geometry

        for i in range(n_edges):
            edge_indexes[:,i] = space.objects_per_dimension[1].object[i].nodes

        for i in range(n_tris):
            tria_indexes[:,i] = space.objects_per_dimension[2].object[i].nodes

        return node_coordinates,edge_indexes,tria_indexes

    def get_xyz(self,node_coordinates,tria_indexes):
        '''
        Parse out triangle corner coordinates. Assume one based indexing in the input (first = 1).

        Returns  xyz[ntri,3,3], where the last dimension [ntri,ncorn,:]  is x,y,z. 


        '''
        ntri =  tria_indexes.shape[1]

        xyz = np.zeros( (ntri, 3, 3 ),dtype=float )

        for itri in range( ntri ):
            for i in range(3):
                xyz[itri,i,:] = node_coordinates[:,tria_indexes[i,itri]-1]

        return xyz
