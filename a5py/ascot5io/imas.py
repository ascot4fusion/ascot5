import numpy as np
import scipy.constants as constants
from a5py.physlib import species, pol2cart, cart2pol_vec
from a5py.ascot5io import options
import unyt
import warnings
from types import SimpleNamespace
import warnings

class a5imas:


    def __init__(self):
        self.ids_name = "ids" #This should be overwritten by dervied classes

    def open(self, user, tokamak, version, shot, run, occurrence=0, backend=None, path=None):

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
                                'backend'    : backend,
                                'path'       : path,
                                'ids_name'   : self.ids_name }

        if hasattr(imas, 'ids'):
            # The "old" 3/4 AL
            self.ids = imas.ids(shot, run)
            self.ids.open_env(user, tokamak, version)

            #ids = ids.split('/')
            #if len(ids) == 1:
            #    occurrence = 0
            #else:
            #    occurrence = int(ids[1])
            #ids = ids[0]


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
        else:
            # The "new" AL 5

            time=0.0

            from imas.imasdef import MDSPLUS_BACKEND
            from imas.imasdef import CLOSEST_SAMPLE
            #self.DB = imas.DBEntry(MDSPLUS_BACKEND, tokamak, shot, run, user_name=user)
            self.DB = imas.DBEntry("imas:"+str(self.ids_coordinates["backend"])
                                   +"?path="+str(self.ids_coordinates["path"]), mode="r")
            self.DB.open()
            self.ids = SimpleNamespace()
            setattr(self.ids,self.ids_name, self.DB.get_slice(self.ids_name, time, CLOSEST_SAMPLE))

        return self.ids

    def setIds(self,ids,ids_coordinates=None):
        '''
        This setter function is intended to be used when the IDS data is
        received as a python object, not read from file.

        ids_coordinates could be e.g. {'user' : "username", 'tokamak' : 'test',
                                       'version': "1.2.3", 'shot': '444',
                                       'run' : "66", 'occurrence' : "0",
                                       'ids_name': "equilibrium" }
        '''
        self.ids = ids
        if ids_coordinates is not None:
            self.ids_coordinates = ids_coordinates
        else:
            self.ids_coordinates = { 'user' : "undefined", 'tokamak' : 'undefined',
                                     'version': "undefined", 'shot': 'undefined',
                                     'run' : "undefined", 'occurrence' : "undefined",
                                     'ids_name': self.ids_name }

    def read(self, user, tokamak, version, shot, run, occurrence=0, backend=None, path=None, **kwargs):
        """
        Open the IDS database and parse the data and then close it.
        """
        self.open( user, tokamak, version, shot, run, occurrence, backend=backend, path=path )
        result = self.parse(**kwargs)
        self.close()

        return result


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


        # get the right imas datastructure
        if self.ids_name == 'wall':
            #ids = data_entry.get(ids_name,occurrence)
            self.ids = imas.wall()
        elif self.ids_name == 'equilibrium':
            self.ids = imas.equilibrium()
        elif self.ids_name == 'distributions':
            self.ids = imas.distributions()
        else:
            raise NotImplementedError("The ids {} has not been implemented in ascot5 imas.py create().".format(self.ids_name))
            # Probably easy enough to implement another ids "newids":
            # "ids = imas.newids()"

    def write_data_entry(self):

        self.data_entry.put(self.ids,self.occurrence)

        self.data_entry.close()


    def close(self):
        import imas
        if hasattr(imas, 'ids'):
            self.ids.close()
        else:
            self.DB.close()

    def fill_mandatory(self,time=[0.0]):
        # mandatory
        self.ids.ids_properties.homogeneous_time = 1
        self.ids.time = np.array(time)


    def description_string(self):
        c = self.ids_coordinates
        s = "DB:{} USER:{} IDS:{} SHOT:{} RUN:{} OCCUR:{}"
        return s.format(c['tokamak'], c['user'], c['ids_name'], c['shot'], c['run'], c['occurrence'])


    def fill_species(self,target_species,anum,znum,charge):

        sp = species.autodetect(anum, znum, charge)

        if sp['name']=='e':
            target_species.type.index         = 1
            target_species.type.name          = "electron"
            target_species.type.description   = "Electron"
        elif round(charge) == 0:
            target_species.type.index         = 4
            target_species.type.name          = "neutral"
            target_species.type.description   = "Neutral species in a single/average state; refer to neutral-structure"
        else:
            target_species.type.index         = 2
            target_species.type.name          = "ion"
            target_species.type.description   = "Ion species in a single/average state; refer to ion-structure"

        if target_species.type.index == 1 or target_species.type.index == 2:
            # ions or electrons

            # We assume single nucleus for our ions. (len(element)==1 && atoms_n==1)
            target_species.ion.element.resize(1)
            target_species.ion.element[0].a       = float(sp['mass'])
            target_species.ion.element[0].z_n     = float(sp['znum'])
            target_species.ion.element[0].atoms_n = 1
            target_species.z_ion                  = float(sp['charge'])
        else:
            # We assume single nucleus for our neutrals. (len(element)==1 && atoms_n==1)
            target_species.neutral.element.resize(1)
            target_species.neutral.element[0].a       = float(sp['mass'])
            target_species.neutral.element[0].z_n     = float(sp['znum'])
            target_species.neutral.element[0].atoms_n = 1
            target_species.z_ion                      = float(sp['charge'])


    def fill_code(self,target_code,runobject):
        """
        Fills in the code version information
        """
        codeversion = runobject.getcodeversion()
        for field in ['name','description','commit','version','repository']:
            setattr( target_code, field, codeversion[field] )

        # List of the code specific parameters in XML format
        # convert parameters into XML
        _, xml_string = options.Opt.schema(runobject.options.read())
        target_code.parameters = xml_string

        #Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used. {dynamic}
        # - output_flag : np.ndarray 1D with int)
        # Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
        target_code.output_flag = np.array( [runobject.get_run_success()] )



    def fill_grid_rz(self,target_grid,r,z):

        target_grid.type.index        = 0
        target_grid.type.name         = "rectangular RZ"
        target_grid.type.description  = "Rectangular grid in the (R,Z) coordinates;"

        target_grid.r = np.array(r) #unyt arrays do not fit into imas
        target_grid.z = np.array(z)

class B_STS(a5imas):
    ''' Read stellarator 3D magnetic field with the conventions laid out in:
        git@github.com:sjjamsa/imas-ggd-b3d.git

        Returns a dict modelled after write_hdf5() in  B_STS.py
    '''

    def __init__(self):
        super().__init__()
        self.ids_name = "equilibrium"



    def parse(self):




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
        self.create( user, tokamak, version, shot, run)
        self.fill(B, metadata)
        self.write_data_entry()

    def fill(self, B, metadata={}):

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


        if 'ids_comment' in metadata :
            self.ids.ids_properties.comment = metadata['ids_comment']
        else:
            self.ids.ids_properties.comment = "equilibrium IDS for testing"

        self.fill_mandatory()


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

        
                                 
class marker(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "distribution_sources"


    def parse(self):
        """
        Read an IMAS 2D-wall into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from
              wall.description_2d[0].limiter.unit[0].outline.r
        and
              wall.description_2d[0].limiter.unit[0].outline.z

        """



        timeIndex = 0
        unit      = 0

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
            out[f] = np.array([],dtype=srcs[0][f].dtype)*srcs[0][f].units
            for s in srcs:
                if s is None:
                    continue
                out[f] = np.concatenate( (out[f],s[f]) )
            try:
                if out[f].units != srcs[0][f].units:
                    print('Warning: unit mismatch')
            except AttributeError:
                # Add the missing units:
                out[f] *= s[f].units

        out['ids']   = np.arange(1,n+1,dtype=int)



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
             units=[unyt.m, unyt.deg,unyt.m,unyt.J,unyt.dimensionless ]
             for i,iind in enumerate(imasinds):
                 out[a5fields[i]] = source.markers[timeIndex].positions[:,np.argwhere(indexes==iind)[0][0]]*units[i]

        elif gyro_type == 1:
            # PARTICLE LOCATION

            # From coordinates:
             #------------------

             # for fortran: use imas_coordinate_identifier, only: IDcoord => coordinate_identifier
             # what about python?
             imasinds = [ 4,    5,    3,  103 ] #,    -1,     -1,  -1 ]  # no vr available!
             a5fields = ['r', 'phi', 'z','vz'] #, 'vr', 'vphi', 'vz']
             units=[unyt.m, unyt.deg,unyt.m,unyt.m/unyt.s ]
             for i,iind in enumerate(imasinds):
                 out[a5fields[i]] = source.markers[timeIndex].positions[:,np.argwhere(indexes==iind)[0][0]]*units[i]
             i_velocity_x, i_velocity_y, i_velocity_z = 101, 102, 103
             vx = source.markers[timeIndex].positions[:,np.argwhere(indexes==i_velocity_x)[0][0]]
             vy = source.markers[timeIndex].positions[:,np.argwhere(indexes==i_velocity_y)[0][0]]
             #vz = source.markers[timeIndex].positions[:,np.argwhere(indexes==i_velocity_z)[0][0]]
             #print('in imas', vx.shape)

             #print("adding vr,vz")
             (x,y,z) = pol2cart(out['r'], out['phi'], out['z'])
             (out['vr'],out['vphi'],_vz) = cart2pol_vec(vx,x,vy,y)
        else:
            raise ValueError("Unknown gyrotype='{}'".format(gyro_type))

        # Weight:
        #--------

        out['weight'] = np.array(source.markers[timeIndex].weights)/unyt.s

        # From species:
        #--------------
        out['anum']  = np.ones_like(out['weight'].v,dtype=int) * int(np.rint(source.species.ion.element[0].a))*unyt.dimensionless
        out['znum']  = np.ones_like(out['weight'].v,dtype=int) * int(np.rint(source.species.ion.element[0].z_n))*unyt.dimensionless

        # Generated:
        #-----------
        out['ids']   = np.arange(1,n+1,dtype=int)*unyt.dimensionless
        out['charge']= np.ones_like(out['weight'].v,dtype=float) * unyt.e* source.species.ion.z_ion
        out['mass']  = np.ones_like(out['weight'].v,dtype=float) * species.autodetect(
            int(source.species.ion.element[0].a),
            int(source.species.ion.element[0].z_n) )["mass"]#/unyt.kg


        # From parameters (outside the source)
        #-------------------------------------
        out['time'] = np.ones_like(out['weight'].v,dtype=int) * time*unyt.s

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

class plasma_1d(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "core_profiles"


    def parse(self,equilibrium_ids=None):
        """
        Read an IMAS 1D-plasma profiles into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from core_profiles.profiles_1d[0]

        The following keys are looked for:
                   nrho, nion, anum, znum, mass, charge, rho,
                   edensity, etemperature, idensity, itemperature

        Ion temperature is read from species-averaged T_i or if that is not available, from the first ion species.

        Rho-tor --> Rho-pol translation is done with values from equilibrium_ids; read that first, e.g.:
           eq=a5py.ascot5io.imas.B_2DS()
           eqdict=eq.read("akaslos","test","3",92436,306)
        """



        timeIndex = 0
        iElement  = 0
        iIonTemp  = 0

        p1d = self.ids.core_profiles.profiles_1d[timeIndex]

        nrho = len(p1d.grid.rho_tor)
        nion = len(p1d.ion)

        anum   = np.zeros(shape=(nion,))
        znum   = np.zeros_like(anum)
        charge = np.zeros_like(anum)
        mass   = np.zeros_like(anum)
        idensity = np.zeros(shape=(nrho,nion))
        ivtor  = np.zeros_like(idensity)
        for i in range(nion):
            znum[i]       = p1d.ion[i].element[iElement].z_n
            anum[i]       = p1d.ion[i].element[iElement].a
            charge[i]     = znum[i] * unyt.e
            #print( species.autodetect( int( anum[i] ), int( znum[i]) ) )
            print( species.autodetect( int( anum[i] ), int( znum[i]) ) )
            mass[i]       = species.autodetect( int( anum[i] ), int( znum[i]) )['mass'] / unyt.amu # Mass should be in AMU
            idensity[:,i] = p1d.ion[i].density_thermal
            if len(p1d.ion[i].rotation_frequency_tor) > 0:
                ivtor[:,i]    = p1d.ion[i].rotation_frequency_tor # [rad.s^-1]

        # Try species-averaged T_i
        if len(p1d.t_i_average) > 0 :
            itemperature = p1d.t_i_average
        else:
            # As a backup solution, use T_i of species 0 (iIonTemp)
            itemperature = p1d.ion[iIonTemp].temperature

        # take the mean of toroidal velocities
        # This is probably all wrong, but as long as they are all zeroes, who cares?
        vtor = ivtor.mean(axis=1)


        edensity     = p1d.electrons.density
        etemperature = p1d.electrons.temperature
        warnings.warn("Reading plasma rotation not yet implemented and is assumed to be zero")
        vtor = edensity * 0

        warnings.warn("Reading plasma rotation not yet implemented and is assumed to be zero")
        vtor = edensity * 0

        rho_tor = p1d.grid.rho_tor

        # https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids2plasmabkg.F90
        # https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids2plasmaEqWallSimu.F90#L281
        if equilibrium_ids is None:
            warnings.warn("Cannot convert rho_pol to rho_tor as no equilibrium ids provided as a parameter")
            p = {
                "nrho"         : nrho,
                "nion"         : nion,
                "anum"         : anum,
                "znum"         : znum,
                "mass"         : mass,
                "charge"       : charge,
                "itemperature" : itemperature,
                "idensity"     : idensity,
                "etemperature" : etemperature,
                "edensity"     : edensity,
                "rho_tor"      : rho_tor,
                "vtor"         : vtor,
            }

            return p


        rho = equilibrium_ids.tor2pol(rho_tor)

        rho = equilibrium_ids.tor2pol(rho_tor)
        p = {
            "nrho"         : nrho,
            "nion"         : nion,
            "anum"         : anum,
            "znum"         : znum,
            "mass"         : mass,
            "charge"       : charge,
            "itemperature" : itemperature,
            "idensity"     : idensity,
            "etemperature" : etemperature,
            "edensity"     : edensity,
            "rho"          : rho,
            "vtor"         : vtor,
        }

        return p


class wall_2d(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "wall"


    def parse(self):
        """
        Read an IMAS 2D-wall into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from
              wall.description_2d[0].limiter.unit[0].outline.r
        and
              wall.description_2d[0].limiter.unit[0].outline.z

        """



        timeIndex = 0
        unit      = 0

        #itm.wall.get()
        r = self.ids.wall.description_2d[timeIndex].limiter.unit[unit].outline.r
        z = self.ids.wall.description_2d[timeIndex].limiter.unit[unit].outline.z


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


    def parse(self):
        """
        Read an IMAS 3D-wall into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from
              wall.description_ggd[0].grid_ggd[0].grid_ggd.space[0]
        Flag is set to the index of description_ggd == 0

        """


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



        # Convert 3D wall data to ASCOT5 format

        xyz = self.get_xyz(node_coordinates,tria_indexes)


        nelements = xyz.shape[0]
        flag      = np.ones( (nelements,1) ) * description_ggd_index
        flagIdList= np.unique(flag)
        flagIdStrings= []
        for f in flagIdList:
            flagIdStrings.append("IMAS 3D-wall ids u:{} db:{} v:{} s:{} r:{} o:{} desc:{} grid:{} sp:{}".format(
                self.ids_coordinates['user'],
                self.ids_coordinates['tokamak'],
                self.ids_coordinates['version'],
                self.ids_coordinates['shot'],
                self.ids_coordinates['run'],
                self.ids_coordinates['occurrence'],
                description_ggd_index,grid_ggd_index,space_index ))
#        print(flagIdList)
#        print(flagIdStrings)
        labels    =  dict( zip(  flagIdStrings, flagIdList ) )


        w = {
            "x1x2x3"         : xyz[:,:,0],
            "y1y2y3"         : xyz[:,:,1],
            "z1z2z3"         : xyz[:,:,2],
            "nelements"      : [[nelements]],
            "flag"           : flag,
            "labels"         : labels,
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
        self.fill(w, metadata)
        self.write_data_entry()

    def fill(self, w, metadata={}):


        n = w['nelements'][0][0]
        xyz=np.zeros( (n, 3, 3) )

        xyz[:,:,0] = w['x1x2x3'][:,:]
        xyz[:,:,1] = w['y1y2y3'][:,:]
        xyz[:,:,2] = w['z1z2z3'][:,:]

        # Save the data
        self.fill_wall_3d_ids(self.ids,xyz,metadata=metadata)


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


class B_2DS(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "equilibrium"


    def parse(self):
        """
        Read an IMAS 2D-equilibrium into a dictionary.
        The dictionary is a drop-in replacement for a dict read from hdf5.

        The following keys are provided:
            rmin, rmax, nr, zmin, zmax, nz,
            axisr, axisz, psi, psi0, psi1,
            (br,bz = will be set to zero,
             bphi contains the full toroidal field)
            psi1 must be exactly the psi at separatrix ( read from global_quantities.psi_boundary )
            psi0 must be exactly the psi at axis, but this can be updated later
            axisr must be exactly the R at axis, but this can be updated later
            axisz must be exactly the z at axis, but this can be updated later

        The data is read from ids.equilibrium

        After parsing, the Object will contain routines to interpolate between rho_pol and rho_tor

        """



        timeIndex = 0
        unit      = 0
        psiscale = 1.0 / ( 2.0 * np.pi )


        # Identify a Rectangular (R,z) grid, index 1
        p2dindex = -1
        for i,profile in enumerate(self.ids.equilibrium.time_slice[timeIndex].profiles_2d):
            if  profile.grid_type.index == 1:
                p2dindex = i
                break
        if p2dindex == -1:
            # let us issue a warning
            import warnings
            warnings.warn("No  Rectangular (R,z) grid found in equilibrium ids profiles_2d")
            return None


        nr   = len(self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].grid.dim1)
        rmin =     self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].grid.dim1[0]
        rmax =     self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].grid.dim1[nr-1]

        nz   = len(self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].grid.dim2)
        zmin =     self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].grid.dim2[0]
        zmax =     self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].grid.dim2[nz-1]

        psi  =     self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].psi      * psiscale
        bphi =     self.ids.equilibrium.time_slice[timeIndex].profiles_2d[p2dindex].b_field_tor

        # These can usually be zero
        br   = np.zeros_like(bphi)
        bz   = np.zeros_like(bphi)

        # values at axis
        axisr=     self.ids.equilibrium.time_slice[timeIndex].global_quantities.magnetic_axis.r
        axisz=     self.ids.equilibrium.time_slice[timeIndex].global_quantities.magnetic_axis.z
        psi0 =     self.ids.equilibrium.time_slice[timeIndex].global_quantities.psi_axis     * psiscale

        # values at separatrix
        psi1  =    self.ids.equilibrium.time_slice[timeIndex].global_quantities.psi_boundary * psiscale

        # for rho_tor -- rho_pol interpolation
        # https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids2plasmaEqWallSimu.F90#L281

        if   len( self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor      ) > 0 :
            rho_tor  = self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor
#        elif len( self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor_norm ) > 0 :
#            import warnings
#            warnings.warn("Using equilibrium.timeslice[].profiles1d.rho_tor_norm (instead of unnormalized)")
#            rho_tor  = self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor_norm
        else:
            rho_tor = None

        psi_prof = self.ids.equilibrium.time_slice[timeIndex].profiles_1d.psi
        self._init_interpolator_functions(rho_tor,psi_prof)



        # for rho_tor -- rho_pol interpolation
        # https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids2plasmaEqWallSimu.F90#L281

        if   len( self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor      ) > 0 :
            rho_tor  = self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor
#        elif len( self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor_norm ) > 0 :
#            import warnings
#            warnings.warn("Using equilibrium.timeslice[].profiles1d.rho_tor_norm (instead of unnormalized)")
#            rho_tor  = self.ids.equilibrium.time_slice[timeIndex].profiles_1d.rho_tor_norm
        else:
            rho_tor = None

        psi_prof = self.ids.equilibrium.time_slice[timeIndex].profiles_1d.psi
        self._init_interpolator_functions(rho_tor,psi_prof)



        b = {
             'rmin' : rmin,
             'rmax' : rmax,
             'nr'   : nr,
             'zmin' : zmin,
             'zmax' : zmax,
             'nz'   : nz,
             'axisr': axisr,
             'axisz': axisz,
             'psi'  : psi,
             'psi0' : psi0,
             'psi1' : psi1,
             'br'   : br,
             'bz'   : bz,
             'bphi' : bphi,
        }

        return b

    def _init_interpolator_functions(self,rho_tor,psi_prof):
        """
        We try to do similar work as this routine:
        https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/rhotorpol.F90
        """


        #TODO, ! If rho_tor doesn't start from 0.000, extrapolate
        if rho_tor[0] > 0:
            warnings.warn("Rho tor-profile does not start from 0.0000, this may be a problem. It was extrapolated in ASCOT4, but has not been implemented in ASCOT5.")


        # Rho_pol as normalized psi:
        rho_pol = np.sqrt(
            (psi_prof     - psi_prof[0] )/
            (psi_prof[-1] - psi_prof[0] ) )

        self.poltor_rho_tor = rho_tor
        self.poltor_rho_pol = rho_pol

        # Add "public" names for the methods as they have been initialized
        self.tor2pol = self._tor2pol
        self.pol2tor = self._pol2tor

    def _tor2pol(self,rho_tor):

        return np.interp( rho_tor, self.poltor_rho_tor, self.poltor_rho_pol)

    def _pol2tor(self,rho_pol):

        return np.interp( rho_pol, self.poltor_rho_pol, self.poltor_rho_tor)


class distributions(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "distributions"

    def write(self, user, tokamak, version, shot, run, runobject, metadata={} ):
        self.create( user, tokamak, version, shot, run)
        self.fill(runobject,metadata)
        self.write_data_entry()

    def fill(self,runobject,metadata={}):
        '''
        The distributions are the '5d' and '5drho distributions you get from the runobject Ascot.data.run_xxxxxxxxx.getdist('5d')
        Thus, e.g. runobject=Ascot.data.run_xxxxxxxxx or a vrun
        '''
        timeIndex = 0

        if 'ids_comment' in metadata :
            self.ids.ids_properties.comment = metadata['ids_comment']
        else:
            self.ids.ids_properties.comment = "distributions IDS for testing"

        self.fill_mandatory()

        self.fill_code(self.ids.code,runobject)

        species = runobject.getspecies()
        anum = float(species['anum'])
        znum = float(species['znum'])

        # Defines how to interpret the spatial coordinates:
        #    1 = given at the actual particle birth point;
        #    2 = given at the gyro centre of the birth point {constant}
        # Corresponds to runobject.getsimmode()=1 --> gyro_type=1
        #                runobject.getsimmode()=2 --> gyro_type=2

        warnings.warn('Not implemented yet [TODO_KONSTA<--SIMPPA]')
        #gyro_type = runobject.getsimmode()
        gyro_type = 1
        if not (gyro_type == 1 or gyro_type == 2) :
            raise ValueError("Unsupported gyro_type from runobject.getsimmode(): '{}' (should be 1 or 2).".format(gyro_type))


        if '5d' in runobject.getdist_list()[0]:
            d5d = runobject.getdist('5d')
            charges5d = d5d.abscissa('charge')
        else:
            charges5d = []


        warnings.warn("SKIPPING rho5d for debugging")
        if '******rho5d******' in runobject.getdist_list()[0]:
            drho5d = runobject.getdist('rho5d')
            chargesrho5d = drho5d.abscissa('charge')
        else:
            chargesrho5d = []


        # To append the distributions one after each other
        distoffset = 0


        ###########
        # dist 5d #
        ###########
        for iDistribution in range(len(charges5d)):

            self.ids.distribution.resize( 1 + iDistribution + distoffset ) # indexing starts from 0, but counting from 1, thus +1
            d = self.ids.distribution[ iDistribution + distoffset ]


            # If is_delta_f=1, then the distribution represents the deviation from a Maxwellian;
            # is_delta_f=0, then the distribution represents all particles, i.e. the full-f solution {constant}
            d.is_delta_f = 1

            d.gyro_type = gyro_type

            charge = charges5d[iDistribution]
            self.fill_species(d.species,anum,znum,charge)

            prof2d = d.profiles_2d.resize(timeIndex+1)
            prof2d = d.profiles_2d[timeIndex]
            self.fill_grid_rz( prof2d.grid, r=d5d.abscissa('r'), z=d5d.abscissa('z') )

            ascot_names= ['density','toroidalcurrent','pressure','electronpowerdep']
            imas_names = ['density_fast']#,'current_fast_phi','pressure_fast'] #,'collisions.electrons.powerthermal']

            # Set each one of the profiles2d
            for iname in range(len(imas_names)):
                moment   = runobject.getdist_moments( runobject.getdist('5d'),ascot_names[iname])
                ordinate = moment.ordinate(ascot_names[iname], toravg=True) # average in phi-direction
                setattr(prof2d,imas_names[iname],np.array(ordinate))


            warningtext = (
                "5D distribution output still WIP; fill in e.g. profiles_2d:\n"
                "density --> density_fast\n"
                "toroidalcurrent --> current_fast_phi\n"
                "pressure --> pressure_fast\n"
                "electronpowerdep --> collisions.electrons.powerthermal\n"
                "\n"
                "The following may not be useful:\n"
                "parallelcurrent : Parallel current\n"
                "chargedensity : Charge density\n"
                "energydensity : Energy density\n"
                "powerdep : Total deposited power\n"
                "ionpowerdep : Power deposited to ions (should be per species)\n"
            )
            warnings.warn(warningtext)
        distoffset += len(charges5d)

        ##############
        # dist rho5d #
        ##############
        warnings.warn("Rho distribution output not implemented")

        for iDistribution in range(len(chargesrho5d)):

            d = self.ids.distribution[ iDistribution + distoffset ]

        distoffset += len(chargesrho5d)
