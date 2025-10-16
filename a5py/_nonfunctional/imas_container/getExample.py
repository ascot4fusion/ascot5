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
