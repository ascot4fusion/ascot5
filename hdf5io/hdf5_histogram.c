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

   Here is an example of the expected order of indexing in ordinate. (The Fortran/matlab order, I guess)
   ordinate = (double*) malloc(abscissaNslots[0] * abscissaNslots[1] * ordinateLength * sizeof(double) );
   for (iDim2 = 0; iDim2 < abscissaNslots[1]; iDim2++ )
     for ( iDim1 = 0; iDim1 < abscissaNslots[0]; iDim1++)
       for ( iOrd = 0; iOrd< ordinateLength; iOrd++)
         ordinate[ iDim2 * ordinateLength * abscissaNslots[0]  +
           iDim1 * ordinateLength +
           iOrd         ]  =  100*(1+iOrd) + 10 * (1+iDim1) +iDim2 +1;

   @param fileName hdf5 file to write to
   @param runpath path to run group histogram being written belongs to
   @param title histogram title
   @param abscissaDim number of dimensions of the abscissa
   @param ordinateDim number of different quantities that are in this histogram
   @param abscissaNslots how many slots in each abscissa dimension (integer array, length abscissaDim)
   @param abscissaMin the lowest abscissa edge for each abscissa dimension (double array, length abscissaDim)
   @param abscissaMax the highest abscissa edge for each abscissa dimension (double array, length abscissaDim)
   @param abscissaUnits a vector of strings (abscissaDim) defining unit for each abscissa dimension
   @param abscissaNames a vector of strings (abscissaDim) defining name for each abscissa dimension
   @param ordinateUnits a vector of strings (ordinateLength) defining unit for each ordinate
   @param ordinateNames vector of strings (ordinateLength) defining name for each ordinate
   @param ordinate a double vector, size = ordinateLength * product(absicssaNslots)

   @return end status, zero on success
*/
int hdf5_histogram_write_uniform_double(
                    const char *fileName,
                    const char *runpath,
                    const char *title,
                    int abscissaDim,
                    int ordinateDim,
                    int *abscissaNslots,
                    double *abscissaMin,
                    double *abscissaMax,
                    char **abscissaUnits,
                    char **abscissaNames,
                    char **ordinateUnits,
                    char **ordinateNames,
                    double *ordinate) {

    hid_t fileHandle;
    hid_t groupHandleHist;
    char path[256];
    char temp_path[256];
    herr_t err;

    hsize_t size[1]; /* Helper variable */


    /* Open file and create "dist" main group and subgroup for this specific distribution. */

    fileHandle = hdf5_open(fileName);
    if( fileHandle < 0 ) return fileHandle;

    strcpy(path, runpath);
    strcat(path, "dists/");

    err = hdf5_find_group(fileHandle, path);
    if(err) {
        hid_t groupHandleDist = H5Gcreate2(fileHandle, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(groupHandleDist < 0) {
            return -1;
        }
        H5Gclose (groupHandleDist);
    }

    strcat(path, title);
    groupHandleHist = H5Gcreate2(fileHandle, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(groupHandleHist < 0) {
        H5Fclose(fileHandle);
        return groupHandleHist;
    }


    /* Add ordinate */

    hsize_t* dims = malloc(sizeof(hsize_t)*(abscissaDim + 1)); /* For storing number of slots in ordinate matrix */

    dims[0] = ordinateDim;
    for (int i=0; i<abscissaDim; i++) {
        dims[i+1] = abscissaNslots[i];
    }

    err = H5LTmake_dataset_double (groupHandleHist , "ordinate", abscissaDim+1, dims, ordinate );
    if(err){
        printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
        H5Gclose(groupHandleHist);
        free(dims);
        H5Fclose(fileHandle);
        return err;
    }


    /* Add ordinate and abscissa dimensions and size of the dimensions */

    size[0] = 1;
    err = H5LTmake_dataset_int(groupHandleHist , "ordinate_ndim", 1, size, &ordinateDim);
    if(err){
        printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
        H5Gclose(groupHandleHist);
        free(dims);
        H5Fclose(fileHandle);
        return err;
    }

    err = H5LTmake_dataset_int(groupHandleHist , "abscissa_ndim", 1, size, &abscissaDim);
    if(err){
        printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
        H5Gclose(groupHandleHist);
        free(dims);
        H5Fclose(fileHandle);
        return err;
    }

    /* Ordinate names and units */
    for (int i=0; i<ordinateDim; i++) {
        sprintf(temp_path, "ordinate_name_%06d",i+1);
        err =  H5LTmake_dataset_string ( groupHandleHist, temp_path, ordinateNames[i] ) ;
        if(err){
            printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
            H5Gclose(groupHandleHist);
            free(dims);
            H5Fclose(fileHandle);
            return err;
        }
        sprintf(temp_path, "ordinate_unit_%06d",i+1);
        err =  H5LTmake_dataset_string ( groupHandleHist, temp_path, ordinateUnits[i] ) ;
        if(err){
            printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
            H5Gclose(groupHandleHist);
            free(dims);
            H5Fclose(fileHandle);
            return err;
        }


    }

    /* Add abscissae */

    for (int i=0; i<abscissaDim; i++) {

        double* abscissa = (double *) malloc( (dims[i+1]+1) * sizeof(double) );
        for(int j=0; j<dims[i+1]+1; j++) {
            abscissa[j] = abscissaMin[i] + j * ((abscissaMax[i] - abscissaMin[i])
                                                /abscissaNslots[i]);
        }

        sprintf(temp_path, "abscissa_nslot_%06d",i+1);
        size[0] = 1;
        err = H5LTmake_dataset_int(groupHandleHist , temp_path, 1, size, &(abscissaNslots[i]));
        if(err){
            printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
            H5Gclose(groupHandleHist);
            free(dims);
            free(abscissa);
            H5Fclose(fileHandle);
            return err;
        }


        sprintf(temp_path, "abscissa_vec_%06d",i+1);
        size[0] = abscissaNslots[i]+1;
        err = H5LTmake_dataset_double (groupHandleHist , temp_path, 1, size, abscissa );
        if(err){
            printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
            H5Gclose(groupHandleHist);
            free(abscissa);
            free(dims);
            H5Fclose(fileHandle);
            return err;
        }
        free(abscissa);

        sprintf(temp_path, "abscissa_name_%06d",i+1);
        err =  H5LTmake_dataset_string ( groupHandleHist, temp_path, abscissaNames[i] );
        if(err){
            printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
            H5Gclose(groupHandleHist);
            free(dims);
            H5Fclose(fileHandle);
            return err;
        }

        sprintf(temp_path, "abscissa_unit_%06d",i+1);
        err =  H5LTmake_dataset_string ( groupHandleHist, temp_path, abscissaUnits[i] );
        if(err){
            printf("%s:%d Problem with make_dataset: %d.\n",__FILE__,__LINE__,err);
            H5Gclose(groupHandleHist);
            free(dims);
            H5Fclose(fileHandle);
            return err;
        }

    }

    H5Gclose(groupHandleHist);
    free(dims);
    H5Fclose(fileHandle);

    return 0;
}
