/**
 * @file hdf5_histogram.c
 * @brief HDF5 Histogram IO
 */
#include "hdf5_helpers.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string.h>
#include <stdlib.h>



/**
   @brief Writes a histogram (ASCOT4 style, version 3) to a hdf5 file. 

   The slot edges would be generated in matlab thus: 
      linspace(abscissaMin(i),abscissaMax(i),abscissaNslots(i)+1)

   @param fileName hdf5 file to write to
   @param title Histogram title
   @param abscissaDim number of dimensions of the abscissa
   @param ordinateLength Number of different quantities that are in this histogram
   @param abscissaNslots how many slots in each abscissa dimension (integer array, length abscissaDim)
   @param abscissaMin the lowest abscissa edge for each abscissa dimension (double array, length abscissaDim)
   @param abscissaMax the highest abscissa edge for each abscissa dimension (double array, length abscissaDim)
   @param abscissaNames a vector of strings (abscissaDim) defining name for each abscissa dimension
   @param abscissaUnits a vector of strings (abscissaDim) defining unit for each abscissa dimension
   @param ordinateNames a vector of strings (ordinateLength) defining name for each ordinate
   @param ordinateUnits a vector of strings (ordinateLength) defining unit for each ordinate
   @param ordinate, a double vector, size = ordinateLength * product(absicssaNslots)
   
   Here is an example of the expected order of indexing in ordinate. (The Fortran/matlab order, I guess)
   ordinate = (double*) malloc(abscissaNslots[0] * abscissaNslots[1] * ordinateLength * sizeof(double) );
   for (iDim2 = 0; iDim2 < abscissaNslots[1]; iDim2++ )
    for ( iDim1 = 0; iDim1 < abscissaNslots[0]; iDim1++)
     for ( iOrd = 0; iOrd< ordinateLength; iOrd++)
      ordinate[ iDim2 * ordinateLength * abscissaNslots[0]  +
                iDim1 * ordinateLength + 
                iOrd         ]  =  100*(1+iOrd) + 10 * (1+iDim1) +iDim2 +1;

   @return end status, zero on success
*/

int hdf5_histogram_write_uniform_double(
		    const char *fileName,
		    const char *runpath, 
		    const char *title,
		    int abscissaDim,
		    int ordinateLength,
		    int *abscissaNslots,
		    double *abscissaMin,
		    double *abscissaMax,
		    char **abscissaUnits,
		    char **abscissaNames,
		    char **ordinateUnits,
		    char **ordinateNames,
		    double *ordinate
			){
  hid_t fileHandle;
  hid_t groupHandleHist,groupHandleOrdinates,groupHandleAbscissae;
  char path[256];
  char temp_path[256];
  herr_t retVal;
  int formatVersion = 3;
  int rank,i,j;
  hsize_t dims[7];
  int slots[6];
  double *abscissa;

  H5open();

  fileHandle = hdf5_open(fileName);
  if( fileHandle < 0 ) return fileHandle;

  strcpy(path, runpath);
  strcat(path, "dists/");

  groupHandleHist = H5Gcreate2(fileHandle, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose (groupHandleHist);
  if(groupHandleHist < 0) {
      return -1;
  }
  
  strcat(path, title);
  groupHandleHist = H5Gcreate2(fileHandle, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //H5Gclose (groupHandleHist);
  if(groupHandleHist < 0) {
    H5Fclose(fileHandle);
    return groupHandleHist;
  }


  // Add ordinate
  rank = 7;
  if( 1 ) {
    dims[6] = ordinateLength;
    for ( i=0;            i<abscissaDim; ++i ) dims[6-(i+1)] = abscissaNslots[i] ;
    //for ( i=abscissaDim;  i<7;           ++i ) dims[6-(i+1)] = 1; This line would just rewrite ordinateLength...
  }else{
    dims[0] = ordinateLength;
    for ( i=0;            i<abscissaDim; ++i ) dims[  (i+1)] = abscissaNslots[i] ;
    for ( i=abscissaDim;  i<7;           ++i ) dims[  (i+1)] = 1;

  }

  retVal = H5LTmake_dataset_double (groupHandleHist , "ordinate", rank, dims, ordinate );
  if(retVal){
    printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,retVal);
    H5Fclose(fileHandle);
    return retVal;
  }


  // Add number of dims, slots, vector length

  retVal = H5LTset_attribute_int(groupHandleHist,"ordinate","number_of_dims", &abscissaDim, 1);
  if(retVal){
    H5Fclose(fileHandle);
    return retVal;
  }

  for ( i=0;            i<abscissaDim; ++i ) slots[i] = abscissaNslots[i] ;
  for ( i=abscissaDim;  i<6;           ++i ) slots[i] = 1;
  retVal = H5LTset_attribute_int(groupHandleHist,"ordinate","number_of_slots", slots, 6);
  if(retVal){
      H5Fclose(fileHandle);
      return retVal;
  }
  retVal = H5LTset_attribute_int(groupHandleHist,"ordinate","vector_length", &ordinateLength, 1);
  if(retVal){
      H5Fclose(fileHandle);
      return retVal;
  }




  // Add ordinate information
  char path_ordinates[256];
  strcpy(path_ordinates, path);
  strcat(path_ordinates,"/ordinates");
  groupHandleOrdinates =  hdf5_create_group( groupHandleHist, "ordinates");
  if(groupHandleOrdinates < 0) {
      H5Fclose(fileHandle);
      return groupHandleOrdinates;
  }
  


  for ( i=0; i<ordinateLength; ++i) {
      sprintf(temp_path, "name_%06d",i+1);
      retVal =  H5LTmake_dataset_string ( groupHandleOrdinates, temp_path, ordinateNames[i] ) ;
      if(retVal){
	  H5Fclose(fileHandle);
	  return retVal;
      }
      sprintf(temp_path, "unit_%06d",i+1);
      retVal =  H5LTmake_dataset_string ( groupHandleOrdinates, temp_path, ordinateUnits[i] ) ;
      if(retVal){
	  H5Fclose(fileHandle);
	  return retVal;
      }
      
      
  }
  H5Gclose(groupHandleOrdinates);



  // Add Abscissae
  char path_abscissae[256];
  strcpy(path_abscissae, path);
  strcat(path_abscissae,"/abscissae");
  groupHandleAbscissae =  hdf5_create_group( groupHandleHist, "abscissae");
  if(groupHandleAbscissae < 0) {
      H5Fclose(fileHandle);
      return groupHandleAbscissae;
  }


  for ( i=0; i<abscissaDim; ++i) {

    rank=1;
    dims[0]=abscissaNslots[i]+1;

    abscissa=(double *) malloc( dims[0] * sizeof(double) );
    for( j=0; j<dims[0]; ++j) abscissa[j] = abscissaMin[i] + j * (( abscissaMax[i] - abscissaMin[i] ) /abscissaNslots[i] );
    sprintf(temp_path, "dim_%06d",i+1);
    retVal = H5LTmake_dataset_double (groupHandleAbscissae , temp_path, rank, dims, abscissa );
    free(abscissa);
    if(retVal){
      printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,retVal);
      H5Fclose(fileHandle);
      return retVal;
    }

    sprintf(temp_path, "name_%06d",i+1);
    retVal =  H5LTmake_dataset_string ( groupHandleAbscissae, temp_path, abscissaNames[i] );
    if(retVal){
      printf("%s:%d Problem with attribute setting: %d.\n",__FILE__,__LINE__,retVal);
      H5Fclose(fileHandle);
      return retVal;
    }

    sprintf(temp_path, "unit_%06d",i+1);
    retVal =  H5LTmake_dataset_string ( groupHandleAbscissae, temp_path, abscissaNames[i] );
    if(retVal){
      printf("%s:%d Problem with attribute setting: %d.\n",__FILE__,__LINE__,retVal);
      H5Fclose(fileHandle);
      return retVal;
    }

  }
  H5Gclose(groupHandleAbscissae);
  H5Gclose(groupHandleHist);



  H5Fclose(fileHandle);
  H5close();

  return 0;
}

