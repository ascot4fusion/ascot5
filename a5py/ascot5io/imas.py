import numpy as np
import scipy.constants as constants
from a5py.physlib import species, pol2cart, cart2pol_vec

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

        self.ids_coordinates = {'user'       : user,
                                'tokamak'    : tokamak,
                                'version'    : version,
                                'shot'       : shot,
                                'run'        : run,
                                'occurrence' : occurrence,
                                'ids_name'   : self.ids_name }

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
        elif self.ids_name == 'equilibrium':
            self.ids = imas.equilibrium()
        else:
            raise NotImplementedError("The ids {} has not been implemented in ascot5 imas.py create().".format(self.ids_name))
            # Probably easy enough to implement another ids "newids":
            # "ids = imas.newids()"

    def write_data_entry(self):
                                       
        self.data_entry.put(self.ids,self.occurrence)

        self.data_entry.close()

    
    def close(self):
        self.ids.close()


    def description_string(self):
        c = self.ids_coordinates
        s = "DB:{} USER:{} IDS:{} SHOT:{} RUN:{} OCCUR:{}"
        return s.format(c['tokamak'], c['user'], c['ids_name'], c['shot'], c['run'], c['occurrence'])


class B_STS(a5imas):
    ''' Read stellarator 3D magnetic field with the conventions laid out in:
        git@github.com:sjjamsa/imas-ggd-b3d.git

        Returns a dict modelled after write_hdf5() in  B_STS.py 
    '''

    def __init__(self):
        super().__init__()
        self.ids_name = "equilibrium"


    def read(self, user, tokamak, version, shot, run, occurrence=0 ):


        # Put this inside the function, not to disturb usage where imas is not available
        import imas

        itm = self.open( user, tokamak, version, shot, run, occurrence )



        time_slice = 0
        ggd_index  = 0
        grids_ggd_index = 0

        psi0 = 0.0 # We assume so
        psi1 = 1.0 # We assume so


        ggd  = self.ids.equilibrium.time_slice[time_slice].ggd[ggd_index]
        grid = self.ids.equilibrium.grids_ggd[grids_ggd_index].grid[0] # Grid for cylindrical R phi z grids. Assumed (by convention) to be in index 0.


        # Check we have R,phi,z grid:
        if grid.identifier.index != 10:
            raise ValueError("Expecting index 10 in grid identifier, instead of {}.".format(grid.identifier.index))
        if len(grid.space) != 3 :
            raise ValueError("Should be 3-dimensional grid.")
        if grid.space[0].coordinates_type[0] != 4 :
            raise ValueError("Coordinate type mismatch [R]")
        if grid.space[1].coordinates_type[0] != 6 :
            raise ValueError("Coordinate type mismatch [phi]")
        if grid.space[2].coordinates_type[0] != 5 :
            raise ValueError("Coordinate type mismatch [z]")


        nR   = len(grid.space[0].objects_per_dimension[0].object)
        nphi = len(grid.space[1].objects_per_dimension[0].object)
        nz   = len(grid.space[2].objects_per_dimension[0].object)
        if nR   < 1 :
            raise ValueError("R is zero length")
        if nphi < 1 :
            raise ValueError("phi is zero length")
        if nz   < 1 :
            raise ValueError("z is zero length")


        R   = np.zeros(shape=(  nR, ) )
        phi = np.zeros(shape=(nphi, ) )
        z   = np.zeros(shape=(  nz, ) )

        for i in range(nR):
            R[i]   = grid.space[0].objects_per_dimension[0].object[i].geometry[0]
        for i in range(nphi):
            phi[i] = grid.space[1].objects_per_dimension[0].object[i].geometry[0]
        for i in range(nz):
            z[i]   = grid.space[2].objects_per_dimension[0].object[i].geometry[0]

        shape = (nR, nphi, nz)
        order = 'F'

        ldata = nR * nphi * nz
        if  len( ggd.b_field_r[   0 ].values) != ldata :
             raise ValueError("B_R size does not match R,phi,z size") 
        if  len( ggd.b_field_tor[   0 ].values) != ldata :
             raise ValueError("B_tor size does not match R,phi,z size") 
        if  len( ggd.b_field_z[   0 ].values) != ldata :
             raise ValueError("B_z size does not match R,phi,z size") 

        B_R   = np.reshape( ggd.b_field_r[   0 ].values, newshape=shape, order=order )
        B_tor = np.reshape( ggd.b_field_tor[ 0 ].values, newshape=shape, order=order )
        B_z   = np.reshape( ggd.b_field_z[   0 ].values, newshape=shape, order=order )

        if len( ggd.psi ) > 0 :
            if  len( ggd.psi[   0 ].values) != ldata :
                raise ValueError("psi data size does not match R,phi,z size") 
            psi_arr   = np.reshape( ggd.psi[   0 ].values, newshape=shape, order=order )
        else:
            psi_arr   = None

        if len( ggd.phi) > 0 :
            if  len( ggd.phi[   0 ].values) != ldata :
                raise ValueError("phi data size does not match R,phi,z size") 
            phi_arr   = np.reshape( ggd.phi[   0 ].values, newshape=shape, order=order )
        else:
            phi_arr   = None

        if len( ggd.theta) > 0 :
            if  len( ggd.theta[   0 ].values) != ldata :
                raise ValueError("theta data size does not match R,phi,z size") 
            theta_arr = np.reshape( ggd.theta[ 0 ].values, newshape=shape, order=order )
        else:
            theta_arr = None

        # Read the magnetix axis from grid[1] (index is "1" by convention)
        if len (   self.ids.equilibrium.grids_ggd[grids_ggd_index].grid    ) > 1 :
            grid = self.ids.equilibrium.grids_ggd[grids_ggd_index].grid[1]
            axis_R = None
            axis_z = None
            nphi = len(grid.space[0].objects_per_dimension[0].object)
            if nphi < 1 :
                raise ValueError("axis data is zero length")

            axis_R   = -1.0 * np.ones( shape=(nphi,) )
            axis_z   = -1.0 * np.ones( shape=(nphi,) )
            axis_phi = -1.0 * np.ones( shape=(nphi,) )

            if (
                    grid.space[0].coordinates_type[0] != 4 or  # R
                    grid.space[0].coordinates_type[1] != 6 or  # phi
                    grid.space[0].coordinates_type[2] != 5 ):  # z
                raise ValueError("expected coordinates_type=[4 6 5]")



            for i in range(nphi):
                axis_R[i]   = grid.space[0].objects_per_dimension[0].object[i].geometry[0]
                axis_phi[i] = grid.space[0].objects_per_dimension[0].object[i].geometry[1]
                axis_z[i]   = grid.space[0].objects_per_dimension[0].object[i].geometry[2]
        else:
            axis_R = None
            axis_z = None


        # Create the dictionary



        B = {}

        Rmin = R[ 0]
        Rmax = R[-1]
        nR   = len(R)
        Pmin = phi[ 0]
        Pmax = phi[-1]
        nP   = len(phi)
        zmin = z[ 0]
        zmax = z[-1]
        nz   = len(z)

        # For the 3D-arrays, the dimension order is (r,phi,z)
        # Check that initial data is the shape we expect it to be
        if (
                B_R.shape[0] != nR or
                B_R.shape[1] != nP or
                B_R.shape[2] != nz ):
            raise ValueError("unexpected array dimensions")





        B["b_rmin"]     =      np.array([Rmin   ])
        B["b_rmax"]     =      np.array([Rmax   ])
        B["b_nr"]       =      np.array([nR     ])
        B["b_phimin"]   =      np.array([Pmin   ])
        B["b_phimax"]   =      np.array([Pmax   ])
        B["b_nphi"]     =      np.array([nP     ])
        B["b_zmin"]     =      np.array([zmin   ])
        B["b_zmax"]     =      np.array([zmax   ])
        B["b_nz"]       =      np.array([nz     ])
        B["psi_rmin"]   =      np.array([Rmin   ])
        B["psi_rmax"]   =      np.array([Rmax   ])
        B["psi_nr"]     =      np.array([nR     ])
        B["psi_phimin"] =      np.array([Pmin   ])
        B["psi_phimax"] =      np.array([Pmax   ])
        B["psi_nphi"]   =      np.array([nP     ])
        B["psi_zmin"]   =      np.array([zmin   ])
        B["psi_zmax"]   =      np.array([zmax   ])
        B["psi_nz"]     =      np.array([nz     ])
        B["axis_phimin"]=      np.array([Pmin   ])
        B["axis_phimax"]=      np.array([Pmax   ])
        B["axis_nphi"]  =      np.array([nP     ])

        B["axisr"]      =      axis_R 
        B["axisz"]      =      axis_z 

        # For the 3D-arrays, the required dimension order is (r,phi,z)
        B["br"]         =      B_R
        B["bphi"]       =      B_tor
        B["bz"]         =      B_z
        B["psi"]        =      psi_arr

        #print('old ', psi_arr.shape)
        #print('new ', B['psi'].shape)
        #print('len ', (nR,nP,nz), 'R,phi,z')

        B["psi0"]       =      np.array([psi0])  # This is a bold assumption atm
        B["psi1"]       =      np.array([psi1])  # This is a bold assumption atm.


        B["desc"]       = self.description_string()

        return B

    def write(self, B, user, tokamak, version, shot, run, metadata={}):

        # Check that all arrays have same size
        if B['b_nr']   != B['psi_nr']:
            raise ValueError('b_nr != psi_nr')
        if B['b_nz']   != B['psi_nz']:
            raise ValueError('b_nz != psi_nz')
        if B['b_nphi'] != B['psi_nphi']:
            raise ValueError('b_nphi != psi_nphi')
        if B['b_nphi'] != B['axis_nphi']:
            raise ValueError('b_nphi != axis_nphi')

        if (
                B['b_rmin'] != B['psi_rmin'] or
                B['b_rmax'] != B['psi_rmax'] or
                B['b_zmin'] != B['psi_zmin'] or
                B['b_zmax'] != B['psi_zmax'] or
                B['b_phimin'] != B['psi_phimin'] or
                B['b_phimax'] != B['psi_phimax'] or
                B['b_phimax'] != B['psi_phimax'] or
                B['b_phimin'] != B['axis_phimin'] or
                B['b_phimax'] != B['axis_phimax']   ):
            raise ValueError('b/psi/axis r/z/phi min/max do not match Expecting identical grids.')
                
        
        index_grids_ggd = 0
        index_time_slice = 0
        nGrids = 2
        index_ggd = 0

        self.create( user, tokamak, version, shot, run)

        if 'ids_comment' in metadata :
            self.ids.ids_properties.comment = metadata['ids_comment']
        else:
            self.ids.ids_properties.comment = "equilibrium IDS for testing"

        # mandatory
        self.ids.ids_properties.homogeneous_time = 1
        self.ids.time = np.array([0.0])

        
        # Fill data

        if len( self.ids.grids_ggd ) <= index_grids_ggd :
            self.ids.grids_ggd.resize(index_grids_ggd+1)
        if len( self.ids.grids_ggd[index_grids_ggd].grid ) < nGrids :
            self.ids.grids_ggd[index_grids_ggd].grid.resize(nGrids)



        # The structured grid
        grid = self.ids.grids_ggd[index_grids_ggd].grid[0]

        if "description_ggd_identifier_name" in metadata :
            grid.identifier.name = metadata[ "description_ggd_identifier_name" ]
        else :
            grid.identifier.name = "structured_spaces"
        if "description_ggd_identifier_desc" in metadata :
            grid.identifier.description = metadata[ "description_ggd_identifier_desc" ]
        else :
            grid.identifier.description = "Cylindrical RphiZ grid"
        grid.identifier.index = 10


        nSpaces = 3 # nodes edges "faces"
        if len(grid.space) < nSpaces:
            grid.space.resize(nSpaces)

        grid.space[0].coordinates_type = np.array( [4] ) # R
        grid.space[1].coordinates_type = np.array( [6] ) # phi
        grid.space[2].coordinates_type = np.array( [5] ) # phi

        objdims=1 # GGD would create additional dimensions, but we don't really need them, do we?
        for iSpace in range(3):
            if len(grid.space[iSpace].objects_per_dimension) < objdims :
                grid.space[   iSpace].objects_per_dimension.resize(objdims)

        nR   = B['b_nr'][0]
        if len(grid.space[0].objects_per_dimension[0].object) != nR :
            grid.space[   0].objects_per_dimension[0].object.resize(nR)
        nPhi = B['b_nphi'][0]
        if len(grid.space[0].objects_per_dimension[0].object) != nPhi :
            grid.space[   1].objects_per_dimension[0].object.resize(nPhi)
        nz   = B['b_nz'][0]
        if len(grid.space[2].objects_per_dimension[0].object) != nz :
            grid.space[   2].objects_per_dimension[0].object.resize(nz)

        R   = np.linspace( B['b_rmin'  ][0], B['b_rmax'  ][0], B['b_nr'  ][0] )
        for i in range(nR):
            grid.space[0].objects_per_dimension[0].object[i].geometry = np.array( [R[i]  ] )
            grid.space[0].objects_per_dimension[0].object[i].nodes = np.array([i+1])
        phi = np.linspace( B['b_phimin'][0], B['b_phimax'][0], B['b_nphi'][0] )
        for i in range(nPhi):
            grid.space[1].objects_per_dimension[0].object[i].geometry = np.array( [phi[i]] )
            grid.space[1].objects_per_dimension[0].object[i].nodes = np.array([i+1])
        z   = np.linspace( B['b_zmin'  ][0], B['b_zmax'  ][0], B['b_nz'  ][0] )
        for i in range(nz):
            grid.space[2].objects_per_dimension[0].object[i].geometry = np.array( [z[i]  ] )
            grid.space[2].objects_per_dimension[0].object[i].nodes = np.array([i+1])
            
        # Todo... edges, sub-spaces

            
        # The magnetic axis grid
        grid = self.ids.grids_ggd[index_grids_ggd].grid[1]

        grid.identifier.index = 1
        grid.identifier.description = 'magnetic axis'

        
        nSpaces = 1 
        if len(grid.space) < nSpaces:
            grid.space.resize(nSpaces)
       
        grid.space[0].identifier.index = 1
        grid.space[0].coordinates_type = np.array([4, 6, 5])

        # nodes & edges
        objdims = 2

        # nodes
        if len(grid.space[0].objects_per_dimension) < objdims :
            grid.space[0].objects_per_dimension.resize(2)
            
        if len(grid.space[0].objects_per_dimension[0].object) != nPhi:
            grid.space[0].objects_per_dimension[0].object.resize(nPhi)

        for i in range(nPhi):
            grid.space[0].objects_per_dimension[0].object[i].geometry = np.array( [  B["axisr"][i], phi[i], B["axisz"][i] ] )
            grid.space[0].objects_per_dimension[0].object[i].nodes = np.array([i+1])

        # edges

        if len(grid.space[0].objects_per_dimension[1].object) != 1:
            grid.space[0].objects_per_dimension[1].object.resize(1)

        grid.space[0].objects_per_dimension[1].object[0].nodes = np.arange( 1, nPhi+1, dtype=int )
            

        if len(grid.grid_subset) < 1:
            grid.grid_subset.resize(1)

        grid.grid_subset[0].identifier.index = 100
        grid.grid_subset[0].dimension = 2

        if len(grid.grid_subset[0].element) < 1:
            grid.grid_subset[0].element.resize(1)
        if len(grid.grid_subset[0].element[0].object) < 1:
            grid.grid_subset[0].element[0].object.resize(1)
    
        grid.grid_subset[0].element[0].object[0].space     = 1
        grid.grid_subset[0].element[0].object[0].dimension = 2
        grid.grid_subset[0].element[0].object[0].index     = 1


        # The fields
        if len( self.ids.time_slice) <= index_time_slice :        
            self.ids.time_slice.resize(index_time_slice+1)
        if len(self.ids.time_slice[index_time_slice].ggd) <= index_ggd :
            self.ids.time_slice[index_time_slice].ggd.resize(index_ggd+1)
        ggd = self.ids.time_slice[index_time_slice].ggd[index_ggd]

        if len(ggd.psi) < 1:
            ggd.psi.resize(1)
        ggd.psi[0].values = B['psi'].flatten(order='F')

        if len(ggd.b_field_r) < 1:
            ggd.b_field_r.resize(1)
        ggd.b_field_r[0].values = B['br'].flatten(order='F')

        if len(ggd.b_field_z) < 1:
            ggd.b_field_z.resize(1)
        ggd.b_field_z[0].values = B['bz'].flatten(order='F')

        if len(ggd.b_field_tor) < 1:
            ggd.b_field_tor.resize(1)
        ggd.b_field_tor[0].values = B['bphi'].flatten(order='F')

        
        self.write_data_entry()
                                 
class marker(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "distribution_sources"


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

        prts = []

        print("Found {} sources.".format(len(self.ids.distribution_sources.source) ) )
        for isrc,source in enumerate(self.ids.distribution_sources.source):
            prts.append(self.read_source(source))


        # Combine prts into a single dict
        out = self.combine_sources(prts)


        return out

    def combine_sources(self, srcs):
        n = 0
        gyro_types = []
        for s in srcs:
            if s is None:
                continue
            n += s['n']
            if 'energy' in s:
                gyro_types.append(1)
            elif 'vr'   in s:
                gyro_types.append(2)
            else:
                gyro_types.append(0)
        print("Total number of markers read: {}".format(n))

        gyro_types=np.array(gyro_types)
        if   np.all( gyro_types == 1):
            gyro_type=1
            fields = ['r', 'phi', 'z', 'energy', 'pitch']
        elif np.all( gyro_types == 2):
            gyro_type=2
            fields = ['r', 'phi', 'z','vz','vr','vphi']
        else:
            raise(ValueError('Unknown gyro type'))

        fields += ['weight','anum','znum','ids','charge','mass','time']

        out={'n':n}
        for f in fields:
            out[f] = np.array([],dtype=srcs[0][f].dtype)
            for s in srcs:
                if s is None:
                    continue
                out[f] = np.concatenate( (out[f],s[f]) )

        out['ids']   = np.arange(1,n+1,dtype=int)

        print(srcs[0]['anum'])
        print(out['anum'])


        return out


    def read_source(self,source,time=0.0):
        """
        Read one source to a dictionary
        """

        # We support only ions
        if not ( source.species.type.index == 2 or
                 source.species.type.index == 3 ):
            return None

        if len(source.markers) <= 0:
            return None

        timeIndex = 0

        out = {}

        # number of markers
        n = len( source.markers[timeIndex].weights )
        if n <= 0:
            return None

        gyro_type = source.gyro_type  #1=actual, 2=gyrocentre

        indexes = np.array([x.index for x in source.markers[timeIndex].coordinate_identifier ])

        if gyro_type == 2:
             # GUIDING CENTRE


             # From coordinates:
             #------------------

             # for fortran: use imas_coordinate_identifier, only: IDcoord => coordinate_identifier
             # what about python?
             imasinds = [ 4,    5,    3,    301,      403,  ]
             a5fields = ['r', 'phi', 'z', 'energy', 'pitch']

             for i,iind in enumerate(imasinds):
                 out[a5fields[i]] = source.markers[timeIndex].positions[:,np.argwhere(indexes==iind)[0][0]]

        elif gyro_type == 1:
            # PARTICLE LOCATION

            # From coordinates:
             #------------------

             # for fortran: use imas_coordinate_identifier, only: IDcoord => coordinate_identifier
             # what about python?
             imasinds = [ 4,    5,    3,  103 ] #,    -1,     -1,  -1 ]  # no vr available!
             a5fields = ['r', 'phi', 'z','vz'] #, 'vr', 'vphi', 'vz']
             for i,iind in enumerate(imasinds):
                 out[a5fields[i]] = source.markers[timeIndex].positions[:,np.argwhere(indexes==iind)[0][0]]
             i_velocity_x, i_velocity_y, i_velocity_z = 101, 102, 103
             vx = source.markers[timeIndex].positions[:,np.argwhere(indexes==i_velocity_x)[0][0]]
             vy = source.markers[timeIndex].positions[:,np.argwhere(indexes==i_velocity_y)[0][0]]
             #vz = source.markers[timeIndex].positions[:,np.argwhere(indexes==i_velocity_z)[0][0]]

             #print("adding vr,vz")
             (x,y,z) = pol2cart(out['r'], out['phi'], out['z'])
             (out['vr'],out['vphi'],_vz) = cart2pol_vec(vx,x,vy,y)
             #out['vr']   = vr
             #out['vphi'] = vphi
             #print("added vr,vz")
        else:
            raise ValueError("Unknown gyrotype='{}'".format(gyro_type))

        # Weight:
        #--------

        out['weight'] = np.array(source.markers[timeIndex].weights)


        # From species:
        #--------------
        out['anum']  = np.ones_like(out['weight'],dtype=int) * int(np.rint(source.species.ion.element[0].a))
        out['znum']  = np.ones_like(out['weight'],dtype=int) * int(np.rint(source.species.ion.element[0].z_n))

        # Generated:
        #-----------
        out['ids']   = np.arange(1,n+1,dtype=int)
        out['charge']= np.ones_like(out['weight'],dtype=float) * constants.elementary_charge * source.species.ion.z_ion
        out['mass']  = np.ones_like(out['weight'],dtype=float) * species.autodetect_species(
            int(source.species.ion.element[0].a),
            int(source.species.ion.element[0].z_n) )['mass']

        # From parameters (outside the source)
        #-------------------------------------
        out['time'] = np.ones_like(out['weight'],dtype=int) * time

        out['n'] = n

        return out

    def display(self):

        timeIndex=0

        #print(self.ids.distribution_sources)
        print("Found {} sources.".format(len(self.ids.distribution_sources.source) ) )
        for isrc,source in enumerate(self.ids.distribution_sources.source):
            print("Source {}:".format(isrc))

            # ion = 2,3 neutral = 4,5, 6=neutron
            print("  Species type {}".format(source.species.type.index) )
            if source.species.type.index == 2 or source.species.type.index == 3:
                if len(source.species.ion.element) > 1:
                    print("  molecule")
                    continue
                print("  a  {} z_ion  {}".format(
                    source.species.ion.element[0].a,
                      source.species.ion.z_ion)    )
            if source.species.type.index == 4 or source.species.type.index == 5:
                if len(source.species.neutral.element) > 1:
                    print("  molecule")
                    continue
                print("  a  {} z_ion  {}".format(
                    source.species.neutral.element[0].a,
                    0.0 ) )
            print("  gyro_type {}".format(source.gyro_type) ) #1=actual, 2=gyro

            if len(source.markers) == 0:
                print("  No markers in source.")
                continue
            indexes = np.array([x.index for x in source.markers[timeIndex].coordinate_identifier ])
            print("  ",end="")
            for ident in source.markers[timeIndex].coordinate_identifier:
                #print("  {}: {}".format(ident.index, ident.name))
                print("{}".format(ident.name),end=', ')
            print("\n  {} markers, summed weight {} particles".format(
                len( source.markers[timeIndex].weights),
                sum( source.markers[timeIndex].weights) ) )
            #print("  Total weight {}".format() ))

            r,phi,z,priv = 4,5,3,-1
            print("  R ranging from {}m to {}m.".format(
                np.amin(source.markers[timeIndex].positions[:,np.argwhere(indexes==r)[0][0]]),
                np.amax(source.markers[timeIndex].positions[:,np.argwhere(indexes==r)[0][0]])
            ))

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


