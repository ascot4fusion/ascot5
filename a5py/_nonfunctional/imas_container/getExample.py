import pickle
import imas
import warnings
import numpy as np

with open('dist5d.pickle', 'rb') as f:
   d5d=pickle.load(f)


distIDS = imas.distributions()
dists   = distIDS.distribution
dists.resize(1)
D=dists[0]



def fill_ggd_d5d( d5d, D, itime=0, itasc=0, irefgrid=0 ):
   '''
   Fill in the general grid discription if the distribution IDS.
   parameters:
     d5d      (input,  from e.g. = vrun.getdist('5d') )
     D        (output, from e.g. dists=imas.distributions().distribution; dists.resize(1) ; D=dists[0] )
     itime    (input, optional(default=0) which time index to fill in the IDS (index here starts with 0)
     itasc    (input, optional(default=0) which time index in the ascot distribution use
     irefgrid (input, optional(default=0) which time index in the ascot distribution use (index here starts with 0, will be +1)

   '''

   reference_implementation_notes='''
The following is a list of needed steps as extracted from manual grid building examples
https://git.iter.org/projects/IMEX/repos/ggd/browse/examples/f90/ids_grid_example1_2dstructured_manual.f90
   2. Set up grid
   - we need one space per abscissa
   - each space can be done separately in a subroutine (gridSetupStruct1dSpace_ex1), and needs
      * coordtype
      * vector of coordinate values
     (To be documented here)
   3. grid subsets are used to define where in the grid the data is actually stored.
      (Assume we have n slots in each dimension.)
       - This is a case of first building up indexes of all points in the grid. That is (n+1)**5 points.
         T
       - Next each neighboring pair of points in connected as a line segment. That is probably another n**5
       - For each face we need to connect four lines.
       - For each cell we need to connect six faces
       - For each tesseract we need to connect eight cells
       - For each penteract we need to connect ten tesseracts

   Luckily, it seems this list does not need to be explictly written.
   See  https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-4.0.0/distributions.html
   # in distribution.ggd.grid.gridsubset.identifier

   # I think we need to request continuation of this list:

   index name
   ----------
   1  nodes         : "All nodes (0D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)."
   2  edges         : "All edges (1D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)"
   43 volumes       : "All volumes (3D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)"
             proposed extension:
   n 4-hypervolumes : "All 4d-hypervolumes (4D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)
   m 5-hypervolumes : "All 5d-hypervolumes (5D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)
   o 6-hypervolumes : "All 6d-hypervolumes (6D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)
   p 7-hypervolumes : "All 7d-hypervolumes (7D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)
   q 8-hypervolumes : "All 8d-hypervolumes (8D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure)




   '''


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
   
   warnings.warn("Species loop not yet implemented")
   iCharge=0
   if len(d5d.abscissa_edges( d5d.abscissae[ 6 ] )) > 2:
      warnings.warn("Only using the index {} in charge dimension.".format(iCharge))


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
      # time index 5 itasc
      # charge should go to the different species
   ]
   
   # We have five spaces, one for each abscissa
   nData = 1
   g.space.resize(len(abscissae))
   for i,a in enumerate(abscissae):

      edges = d5d.abscissa_edges( d5d.abscissae[ a[0] ] )

      S = g.space[i]
      
      S.identifier.name         = 'primary_standard'
      S.identifier.index        = 1
      S.identifier.description  = 'Primary space defining the standard grid'

      S.geometry_type.name        = 'standard'
      S.geometry_type.index       = 0
      S.geometry_type.description = 'standard'

      S.coordinates_type.resize(1)
      C = S.coordinates_type[0]
      C.name        = a[1]
      C.index       = a[2]
      C.description = a[3]
  
      S.objects_per_dimension.resize(1)
      S.objects_per_dimension[0].object.resize( len(edges) )
      for j,e in enumerate(edges):
         S.objects_per_dimension[0].object[j].geometry.resize(1)
         S.objects_per_dimension[0].object[j].geometry[0] = e

      nData *= ( len(edges) -1 ) # The data is on the faces


   # Write the data
   A = D.ggd[itime].amplitude

   A.grid_index = irefgrid

      
   A.values.resize(nData)

   # The data organization, copied from https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-4.0.0/distributions.html
   # in distribution.ggd.grid.gridsubset.identifier for index 1: "nodes"

   '''
   All nodes (0D) belonging to the associated spaces, implicit declaration (no need to replicate the grid elements in the grid_subset structure).
   In case of a structured grid represented with multiple 1D spaces,
   the order of the implicit elements in the grid_subset follows Fortran ordering, i.e.
   iterate always on nodes of the first space first,
   then move to the second node of the second space, ... :
   [((s1_1 to s1_end), s2_1, s3_1 ... sN_1), (((s1_1 to s1_end), s2_2, s3_1, ... sN_1)), ... ((s1_1 to s1_end), s2_end, s3_end ... sN_end)]
   '''



   # Here is a leap of fate that the data gets correctly organized:
   A.values[:] = d5d.distribution().value[:,:,:,:,:,itasc,iCharge].ravel()




print('Starting to fill distribution IDS ggd part')

fill_ggd_d5d( d5d, D)

print('done filling ggd')

print('Preparing to write IDS')
'''
   https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation//ACCESS-LAYER-doc/python/5.4/examples/001_3_open_URI_path.html
'''




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
