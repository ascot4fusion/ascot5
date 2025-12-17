import pickle
import imas
import warnings
import numpy as np
import copy
import unyt

with open('dist5d.pickle', 'rb') as f:
   P=pickle.load(f)
d5d  = P['d5d']
mass = P['mass']


distIDS = imas.distributions()
dists   = distIDS.distribution
dists.resize(1)
D=dists[0]



def fill_ggd_d5d( d5d, D, mass, itime=0, irefgrid=0 ):
   """
   Fill in the general grid discription of the distribution IDS.
   parameters:
     d5d      (input,  from e.g. = vrun.getdist('5d') )
     D        (output, from e.g. dists=imas.distributions().distribution; dists.resize(1) ; D=dists[0] )
     itime    (input, optional(default=0) which time index to fill in the IDS (index here starts with 0)
     irefgrid (input, optional(default=0) which time index in the ascot distribution use (index here starts with 0, will be +1)

   """

   # One can compare to the ASCOT4 implementation:
   # The call is                here: https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids2distribution2ids.F90?ref_type=heads#L293
   # But the GGD part is really here: https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids_grid_structured.f90



   # Sanity check:
   if (d5d.abscissae[0] != 'r'     or
       d5d.abscissae[1] != 'phi'   or
       d5d.abscissae[2] != 'z'     or
       d5d.abscissae[3] != 'ppar'  or
       d5d.abscissae[4] != 'pperp' or
       d5d.abscissae[5] != 'time'  or
       d5d.abscissae[6] != 'charge'   ):
      warnings.warn("The 5d-distribution abscissae are not what is expected. Skipping ggd")
      return


   # Let us integrate away charge and time.
   # First make a deepcopy, because the distribution is changed by the call.
   d5d = copy.deepcopy(d5d)
   d5d.integrate(charge=np.s_[:], time=np.s_[:])


   # Make sure we have an allocated ggd
   if len(D.ggd) <= 0:
      D.ggd.resize(1)

   g = D.ggd[itime].grid


   g.identifier.name = "structured_spaces"
   g.identifier.index = 10
   g.identifier.description = "ASCOT5 5D distribution in"
   for dim in d5d.abscissae:
      g.identifier.description += ' {}'.format(dim)

   # Prepare the spaces...abscissaes

   abscissae=[
      # index ascot, name, index, description
      ( 0, 'r',                          4, 'Major radius'),
      ( 1, 'phi',                        5, 'Toroidal angle'),
      ( 2, 'z',                          3, 'Vertical coordinate z'),
      ( 3, 'momentum_parallel',        201, 'Component of the relativistic momentum vector parallel to the magnetic field'),
      ( 4, 'momentum_perpendicular',   202, 'Component of the relativistic momentum vector perpendicular to the magnetic field'),
      # time index 5 integrated away
      # charge       integrated away
   ]

   # We have five spaces, one for each abscissa
   nData = 1
   g.space.resize(len(abscissae))
   for i,a in enumerate(abscissae):

      edges = d5d.abscissa_edges( d5d.abscissae[ a[0] ] )
      grid_setup_struct1d_space( g.space[i], a, edges )

      nData *= ( len(edges) -1 ) # The data is on the faces

   # Write the data
   A = D.ggd[itime].amplitude

   A.grid_index = irefgrid

   # Comparing to the ASCOT4 4D-distribution.
   # https://version.aalto.fi/gitlab/ascot4/ascot4-trunk/-/blob/master/ids/ids_grid_structured.f90?ref_type=heads#L35
   #    vs
   # https://github.com/iterorganization/IMAS-Data-Dictionary/blob/develop/schemas/utilities/ggd_subset_identifier.xml#L70
   #
   # warnings.warn('ASCOT4 uses grid_subset_index "4" for 4d-distribution, and calls it "all 4d objects", while the DD has 200 for 4d. For 5D ascot5 uses 201.')
   # It actually seems that ASCOT4 does not write the subset at all. Nothing seems to set the grid_subset_index.


   A.grid_subset_index = 201 # "all 5-hypervolumes" https://imas-data-dictionary.readthedocs.io/en/latest/generated/identifier/ggd_subset_identifier.html#identifier-utilities-ggd_subset_identifier.xml
                             # See also https://github.com/iterorganization/IMAS-Data-Dictionary/issues/166

   # A.grid_subset_index *= -1 # We are not following the units of the data dictionary

   A.values.resize(nData)

   # The data organization, copied from https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-4.0.0/distributions.html
   # in distribution.ggd.grid.gridsubset.identifier for index 1: "nodes"
   # Or https://imas-data-dictionary.readthedocs.io/en/latest/generated/identifier/ggd_subset_identifier.html#identifier-utilities-ggd_subset_identifier.xml

   """
   All nodes (0D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure).
   In case of a structured grid represented with multiple 1D spaces,
   the order of the implicit elements in the grid_subset follows Fortran ordering, i.e.
   iterate always on nodes of the first space first,
   then move to the second node of the second space, ... :
   [((s1_1 to s1_end), s2_1, s3_1 ... sN_1), (((s1_1 to s1_end), s2_2, s3_1, ... sN_1)), ... ((s1_1 to s1_end), s2_end, s3_end ... sN_end)]
   """


   """
   from the documentation of ravel():
   order {‘C’,’F’, ‘A’, ‘K’}, optional

   The elements of a are read using this index order.

   - ‘C’ means to index the elements in row-major, C-style order, with the last axis index changing fastest, back to the first axis index changing slowest.
   - ‘F’ means to index the elements in column-major, Fortran-style order, with the first index changing fastest, and the last index changing slowest.
   - ‘A’ means to read the elements in Fortran-like index order if a is Fortran contiguous in memory, C-like order otherwise.
   - ‘K’ means to read the elements in the order they occur in memory, except for reversing the data when strides are negative.

   Note that the ‘C’ and ‘F’ options take no account of the memory layout of the underlying array, and only refer to the order of axis indexing.
   By default, ‘C’ index order is used.
   """

   ravel_order = 'F'
   warnings.warn("Assuming that ravel from ASCOT5 5d-distribution to Fortran works with ravel order '{}'. Not checked.".format(ravel_order))

   # Prepare to move from density to number particles per slot. This is the density of each slot:
   V = d5d.phasespacevolume()

   VDD = calculate_phase_space_volume_ppar_pperp(mass=mass,
                                                 r    = d5d.abscissa_edges( 'r'     ),
                                                 phi  = d5d.abscissa_edges( 'phi'   ),
                                                 z    = d5d.abscissa_edges( 'z'     ),
                                                 ppar = d5d.abscissa_edges( 'ppar'  ),
                                                 pperp= d5d.abscissa_edges( 'pperp' )  )

   # The unit is supposed to be m^-6.s^3 ( https://github.com/iterorganization/IMAS-Data-Dictionary/issues/168 )
   A.values[:] = d5d.distribution().value[:,:,:,:,:].ravel(order=ravel_order) * V.ravel(order=ravel_order) / VDD.ravel(order=ravel_order)


#   warnings.warn("Storing the used phase space volume to ggd[itime].grid.grid_subset[0].metric.jacobian (nonstandard usage)")
#   # except we cannot simply do that:
#   # ValidationError: Element 'distribution[0]/ggd[0]/grid/grid_subset[0]/metric/jacobian' must have its coordinate in dimension 1 (any of 'distribution(i1)/ggd(itime)/grid/grid_subset(i2)/element') filled.
#   D.ggd[itime].grid.grid_subset.resize(1)
#   D.ggd[itime].grid.grid_subset[0].metric.jacobian = VDD.ravel(order=ravel_order)

def calculate_phase_space_volume_ppar_pperp(mass,r,phi,z,ppar,pperp):
   """Calculate the phase space volume of each histogram slot in an ascot5 5d distribution of r,phi,z,ppar,pperp

   The volume in m**3 * (m/s)**3

    Parameters
    ----------

"""
   #from a5py.physlib import velocity_momentum
   #print(mass)

   V = np.ones( ( len(r)-1, len(phi)-1, len(z)-1, len(ppar)-1, len(pperp)-1 ) ) *unyt.m*unyt.m*unyt.m * unyt.m*unyt.m*unyt.m /unyt.s/unyt.s/unyt.s

   # Convert momentum into velocity
   vpar = np.zeros_like(ppar.value) * unyt.m/unyt.s
   for i,p in enumerate(ppar):
      vpar[i] = p/mass #velocity_momentum(mass, np.array([p,]))

   # Convert momentum into velocity
   vperp = np.zeros_like(pperp.value) * unyt.m/unyt.s
   for i,p in enumerate(pperp):
      vperp[i] = p/mass #velocity_momentum(mass, p)

   # Calculate volumes in spatial space

   # Radial slice, surface area of a disc, minus surface area of a smaller disc: pi*R_o**2 - pi*R_i**2
   for i in range( len(r)-1 ):
      V[i,:,:,:,:] *= np.pi * ( r[i+1]*r[i+1] - r[i]*r[i] ) / ( unyt.m * unyt.m)

   # toroidal segment
   for i in range( len(phi)-1 ):
      V[:,i,:,:,:] *=  ( phi[i+1] - phi[i] ) /  ( 2 * np.pi * unyt.rad )

   # vertical slices
   for i in range( len(z)-1 ):
      V[:,:,i,:,:] *=  ( z[i+1] - z[i] )  / ( unyt.m )

   # Calculate volumes in velocity space
   # The geometry is analogous to R,phi,z <--> vperp,gyrophase,vpar

   # parallel slices
   for i in range( len(ppar)-1 ):
      V[:,:,:,i,:] *=  ( vpar[i+1] - vpar[i] )  / ( unyt.m / unyt.s )

   # Radial slice, surface area of a disc, minus surface area of a smaller disc: pi*R_o**2 - pi*R_i**2
   for i in range( len(pperp)-1 ):
      V[:,:,:,:,i] *= np.pi * ( vperp[i+1]*vperp[i+1] - vperp[i]*vperp[i] )  / ( unyt.m / unyt.s )**2

   # (We have all gyrophases included, so no gyrophase related operation needed)

   #print('V[0]',V[0,0,0,0,0])

   return V


def grid_setup_struct1d_space(space, coordtype, nodes, periodic=False):
   """
input:
  space = ggd[itime].grid.space[i]
  coordtype = ('ignored', 'name', 'index', 'description' )
  nodes = list of edge coordinates
"""

   print("  Setting space for {}.".format(coordtype[3]) )
   GRID_UNDEFINED = 0  # https://git.iter.org/projects/IMEX/repos/ggd/browse/src/f90/ids_grid_common.f90#26

   S = space

   S.identifier.name         = 'primary_standard'
   S.identifier.index        = 1
   S.identifier.description  = 'Primary space defining the standard grid'

   S.geometry_type.name        = 'standard'
   S.geometry_type.index       = 0
   S.geometry_type.description = 'standard'

   # https://git.iter.org/projects/IMEX/repos/ggd/browse/examples/f90/ids_grid_example1_2dstructured_manual.f90#280

   S.coordinates_type.resize(1)
   C = S.coordinates_type[0]
   C.name        = coordtype[1]
   C.index       = coordtype[2]
   C.description = coordtype[3]

   # We have nodes and edges
   S.objects_per_dimension.resize(1+1)

   # Nodes
   S.objects_per_dimension[0].object.resize( len(nodes) )
   for j,e in enumerate(nodes):
      S.objects_per_dimension[0].object[j].geometry.resize(1)
      S.objects_per_dimension[0].object[j].geometry[0] = e
      S.objects_per_dimension[0].object[j].nodes.resize(1)
      S.objects_per_dimension[0].object[j].nodes[0] = j+1 # fortran indexing starts with 1

   # edges
   # https://git.iter.org/projects/IMEX/repos/ggd/browse/examples/f90/ids_grid_example1_2dstructured_manual.f90#309

   nobjects = len(nodes) - 1
   S.objects_per_dimension[1].object.resize( nobjects )
   for j in range(nobjects):
      S.objects_per_dimension[1].object[j].nodes.resize(2)
      S.objects_per_dimension[1].object[j].boundary.resize(2)
      for i in range(2):
         S.objects_per_dimension[1].object[j].nodes[i] = j+i+1 # fortran indexing starts with 1
         S.objects_per_dimension[1].object[j].boundary[i].index = i+1 # fortran indexing starts with 1
         S.objects_per_dimension[1].object[j].boundary[i].neighbours.resize(1)

      if periodic and j+i+1 > len(nodes) :  # we use the last value i=1
         S.objects_per_dimension[1].object[j].nodes[i] = 1
         S.objects_per_dimension[1].object[j].boundary[i].index = 1

      # edges
      # left neighbour
      S.objects_per_dimension[1].object[j].boundary[0].neighbours[0] = j - 1 + 1
      # right neighbour
      S.objects_per_dimension[1].object[j].boundary[1].neighbours[0] = j + 1 + 1



   # First edge left neighbour gets special treatment
   # https://git.iter.org/projects/IMEX/repos/ggd/browse/examples/f90/ids_grid_example1_2dstructured_manual.f90#379
   if not periodic:
      S.objects_per_dimension[1].object[0].boundary[0].neighbours[0] = GRID_UNDEFINED
   else:
      S.objects_per_dimension[1].object[0].boundary[0].neighbours[0] = nobjects


print('Starting to fill distribution IDS ggd part')

fill_ggd_d5d( d5d, D, mass)

print('done filling ggd')

print('Preparing to write IDS')
"""
   https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation//ACCESS-LAYER-doc/python/5.4/examples/001_3_open_URI_path.html
"""




from a5py.templates.imasinterface import write_ids

# set mandatory field
distIDS.ids_properties.homogeneous_time = imas.imasdef.IDS_TIME_MODE_HOMOGENEOUS

distIDS.time = np.array([0.0])



USE_A5_WRITE = True

backend="mdsplus"
backend="hdf5"

print(f'writing IDS with "{backend}" backend')


#write_ids(ids, query, backend="hdf5", occurrence=0, host=None):
if USE_A5_WRITE:
   write_ids(ids=distIDS, query=f"/ids/testdb_{backend}", backend=backend)
else:
   with imas.DBEntry(f'imas:{backend}?path=/ids/testdb_{backend}', mode="w") as entry:
      entry.put(distIDS)

print('writing complete')
