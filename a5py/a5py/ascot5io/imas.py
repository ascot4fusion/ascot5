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

    def create(self, user, tokamak, version, shot, run, occurrence=0):

        # Put this inside the function, not to disturb usage where imas is not available
        import imas

        
        backend = imas.imasdef.MDSPLUS_BACKEND

        self.occurrence = 0
        self.data_entry = imas.DBEntry(backend_id=backend,
                                  db_name=tokamak, shot=shot, run=run,
                                  user_name=user, data_version=version )
        

        #data_entry.open()
        self.data_entry.create()
                                  
        if self.ids_name == 'wall':
            #ids = data_entry.get(ids_name,occurrence)
            self.ids = imas.wall()
        else:
            raise NotImplementedError("The ids {} has not been implemented in ascot5 imas.py create().".format(self.ids_name))
            # Probably easy enough to implement another ids "newids":
            # "ids = imas.newids()"

    def write_data_entry(self):
                                       
        self.data_entry.put(self.ids,self.occurrence)

        self.data_entry.close()

    
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
            "nelements"      : [[nelements]],
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

    def fill_wall_3d_ids(self,xyz,metadata={}):
        
        # Let's first do the conversion
        nodes,edges,faces = self.xyz_to_nodes_edges_faces(xyz)

        ids = self.ids
        
        if 'ids_comment' in metadata :
            ids.ids_properties.comment = metadata['ids_comment']
        else:
            ids.ids_properties.comment = "wall IDS saved by ASCOT5 imas.py"

        # mandatory
        ids.ids_properties.homogeneous_time = 1
        ids.time = np.array([0.0])


        index_description_ggd = 0
        index_grid_ggd        = 0
        index_space           = 0

        if len( ids.description_ggd ) <= index_description_ggd :
            ids.description_ggd.resize( index_description_ggd + 1 )
        if len( ids.description_ggd[index_description_ggd].grid_ggd ) <= index_grid_ggd :
            ids.description_ggd[index_description_ggd].grid_ggd.resize( index_grid_ggd + 1 )
        if len( ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].space ) <= index_space :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].space.resize(index_space + 1)

        if "description_ggd_identifier_name" in metadata :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.name = metadata[ "description_ggd_identifier_name" ]
        else :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.name = "ascot5wall"
        if "description_ggd_identifier_desc" in metadata :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.description = metadata[ "description_ggd_identifier_desc" ]
        else :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.description = "Wall stored by ASCOT5 imas.py"
        ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.index = 1


        space = ids.\
            description_ggd[index_description_ggd].\
            grid_ggd[index_grid_ggd].\
            space[index_space]

        # x,y,z space
        space.coordinates_type = np.array([1,2,3])

        space.objects_per_dimension.resize(3)

        space.objects_per_dimension[0].object.resize(len(nodes))
        for i,node in enumerate(nodes):
            space.objects_per_dimension[0].object[i].geometry = node
            space.objects_per_dimension[0].object[i].nodes = np.ones( (1,),dtype=int )*(i+1)

        space.objects_per_dimension[1].object.resize(len(edges))
        for i,edge in enumerate(edges):
            space.objects_per_dimension[1].object[i].nodes = edge

        space.objects_per_dimension[2].object.resize(len(faces))
        for i,face in enumerate(faces):
            space.objects_per_dimension[2].object[i].nodes = face

    def xyz_to_nodes_edges_faces(self,xyz):
        '''
        Input must be xyz[ntri,3,3], where the last dimension [ntri,ncorn,:]  is x,y,z. 
        '''

        #remove duplicate nodes?
        remove_duplicate_nodes = True

        remove_duplicate_edges = True
        if remove_duplicate_edges:
            remove_duplicate_nodes = True

        # number of triangles
        n = xyz.shape[0]

        # Every node is independent. There may be duplicates.
        nodes = np.copy(xyz).reshape( (3*n,3) )

        # We have one index per corner = 3*n
        faces = np.arange(0,3*n,dtype=int).reshape(n,3)

        # a row of edges
        a = np.array([ [0,1],[1,2],[2,0] ]).reshape( (1,3,2) )
        # a column of triangle indexes
        b = np.arange(0,3*n,3,dtype=int).reshape((n,1,1))
        edges = a+b
        edges = edges.reshape(3*n,2)

        if remove_duplicate_nodes :
            nodes,reverse =  np.unique(nodes, return_inverse=True,return_index=False,
                                            axis=0)

            edges = reverse[edges]
            faces = reverse[faces]

        if remove_duplicate_edges :
            edges = np.unique(edges,axis=0)

        # remember to edges++ and faces++ to get one-based indexing 
        return nodes+1,edges+1,faces+1

    def write(self, w, user, tokamak, version, shot, run, metadata={}):

        self.create( user, tokamak, version, shot, run)

        n = w['nelements'][0][0]
        xyz=np.zeros( (n, 3, 3) )

        xyz[:,:,0] = w['x1x2x3'][:,:]
        xyz[:,:,1] = w['y1y2y3'][:,:]
        xyz[:,:,2] = w['z1z2z3'][:,:]
                                       
        # Save the data
        self.fill_wall_3d_ids(self.ids,xyz,metadata=metadata)

                                       
        self.write_data_entry()
                                 

    def fill_wall_3d_ids(self,ids,xyz,metadata={}):

        # Let's first do the conversion
        nodes,edges,faces = self.xyz_to_nodes_edges_faces(xyz)

        if 'ids_comment' in metadata :
            ids.ids_properties.comment = metadata['ids_comment']
        else:
            ids.ids_properties.comment = "wall IDS for testing"

        # mandatory
        ids.ids_properties.homogeneous_time = 1
        ids.time = np.array([0.0])


        index_description_ggd = 0
        index_grid_ggd        = 0
        index_space           = 0

        if len( ids.description_ggd ) <= index_description_ggd :
            ids.description_ggd.resize( index_description_ggd + 1 )
        if len( ids.description_ggd[index_description_ggd].grid_ggd ) <= index_grid_ggd :
            ids.description_ggd[index_description_ggd].grid_ggd.resize( index_grid_ggd + 1 )
        if len( ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].space ) <= index_space :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].space.resize(index_space + 1)

        if "description_ggd_identifier_name" in metadata :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.name = metadata[ "description_ggd_identifier_name" ]
        else :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.name = "MyMesh"
        if "description_ggd_identifier_desc" in metadata :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.description = metadata[ "description_ggd_identifier_desc" ]
        else :
            ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.description = "SMITER mesh MyMesh"
        ids.description_ggd[index_description_ggd].grid_ggd[index_grid_ggd].identifier.index = 1


        space = ids.\
            description_ggd[index_description_ggd].\
            grid_ggd[index_grid_ggd].\
            space[index_space]

        # x,y,z space
        space.coordinates_type = np.array([1,2,3])

        space.objects_per_dimension.resize(3)

        space.objects_per_dimension[0].object.resize(len(nodes))
        for i,node in enumerate(nodes):
            space.objects_per_dimension[0].object[i].geometry = node
            space.objects_per_dimension[0].object[i].nodes = np.ones( (1,),dtype=int )*(i+1)

        space.objects_per_dimension[1].object.resize(len(edges))
        for i,edge in enumerate(edges):
            space.objects_per_dimension[1].object[i].nodes = edge

        space.objects_per_dimension[2].object.resize(len(faces))
        for i,face in enumerate(faces):
            space.objects_per_dimension[2].object[i].nodes = face

    def xyz_to_nodes_edges_faces(self,xyz):
        '''
        Input must be xyz[ntri,3,3], where the last dimension [ntri,ncorn,:]  is x,y,z. 
        '''

        #remove duplicate nodes?
        remove_duplicate_nodes = True

        remove_duplicate_edges = True
        if remove_duplicate_edges:
            remove_duplicate_nodes = True

        # number of trianglesgg59
        n = xyz.shape[0]

        # Every node is independent. There may be duplicates.
        nodes = np.copy(xyz).reshape( (3*n,3) )

        # We have one index per corner = 3*n
        faces = np.arange(0,3*n,dtype=int).reshape(n,3)

        # a row of edges
        a = np.array([ [0,1],[1,2],[2,0] ]).reshape( (1,3,2) )
        # a column of triangle indexes
        b = np.arange(0,3*n,3,dtype=int).reshape((n,1,1))
        edges = a+b
        edges = edges.reshape(3*n,2)

        if remove_duplicate_nodes :
            nodes,reverse =  np.unique(nodes, return_inverse=True,return_index=False,
                                            axis=0)

            edges = reverse[edges]
            faces = reverse[faces]

        if remove_duplicate_edges :
            edges = np.unique(edges,axis=0)

        # remember to edges++ and faces++ to get one-based indexing 
        return nodes+1,edges+1,faces+1


