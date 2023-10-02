/**
 * @file hdf5_histogram.c
 * @brief HDF5 Histogram IO
 */
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_histogram.h"

/**
 * @brief Write a histogram with uniform grid to HDF5 file
 *
 * Histogram can have up to 99 dimensions (abscissae) and likewise the ordinate
 * can be multi-dimensional.
 *
 * The ordinate is stored as Ndim array where format is
 * ordinateDim x abscissa01Dim x abscissa02Dim ...
 *
 * The datasets that are written are the ordinate, ordinate dimension, abscissa
 * dimension, abscissa grid vectors (that define edges of each cell), and
 * abscissa grid vector lengths. Abscissae names and units and ordinate names
 * and units are stored as attributes in ordinate and abscissa grid vectors.
 *
 * The ordinate data is expected to be array with format
 *
 * @param f HDF5 file id
 * @param path full path including group name where histogram is stored
 * @param abscissaDim number of abscissa dimensions
 * @param ordinateDim number of ordinate dimensions
 * @param abscissaNslots array with abscissa dimensions
 * @param abscissaMin the lowest abscissa edge for each abscissa dimension
 * @param abscissaMax the highest abscissa edge for each abscissa dimension
 * @param abscissaUnits array with abscissa units
 * @param abscissaNames array with abscissa names
 * @param ordinateUnits array with ordinate units
 * @param ordinateNames array with ordinate names
 * @param ordinate ordinate data
 *
 * @return zero on success
 */
int hdf5_histogram_write_uniform_double(
    hid_t f, const char *path, int abscissaDim, int ordinateDim,
    int *abscissaNslots, double *abscissaMin, double *abscissaMax,
    char **abscissaUnits, char **abscissaNames,
    char **ordinateUnits, char **ordinateNames, double *ordinate) {

    char temppath[256]; /* Helper variable */

    /* Create histogram group */
    hid_t histogram = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(histogram < 0) {
        return 1;
    }

    /* Dimensions of ordinate data to be written */
    hsize_t dims[100];
    dims[0] = ordinateDim;
    for (int i=0; i<abscissaDim; i++) {
        dims[i+1] = abscissaNslots[i];
    }

    /* Write ordinate including its names and units */
    herr_t err;
    err = H5LTmake_dataset_double(histogram , "ordinate", abscissaDim+1,
                                  dims, ordinate);
    if(err){
        H5Gclose(histogram);
        return err;
    }
    for(int i=0; i<ordinateDim; i++) {
        sprintf(temppath, "name_%02d", i);
        H5LTset_attribute_string(histogram, "ordinate", temppath,
                                 ordinateNames[i]);
        sprintf(temppath, "unit_%02d", i);
        H5LTset_attribute_string(histogram, "ordinate", temppath,
                                 ordinateUnits[i]);
    }


    /* Write ordinate and abscissa dimensions */
    hsize_t dimsize = 1;
    err = H5LTmake_dataset_int(histogram , "ordinate_ndim", 1, &dimsize,
                               &ordinateDim);
    if(err){
        H5Gclose(histogram);
        return 1;
    }
    err = H5LTmake_dataset_int(histogram , "abscissa_ndim", 1, &dimsize,
                               &abscissaDim);
    if(err){
        H5Gclose(histogram);
        return 1;
    }

    /* Write abscissae */
    for (int i=0; i<abscissaDim; i++) {

        double* abscissavec = (double *) malloc( (dims[i+1]+1)*sizeof(double) );
        for(int j=0; j<dims[i+1]+1; j++) {
            abscissavec[j] =
                abscissaMin[i] + j * ( (abscissaMax[i] - abscissaMin[i])
                                       / abscissaNslots[i] );
        }

        /* Vector */
        char abscissapath[256];
        sprintf(abscissapath, "abscissa_vec_%02d",i+1);
        hsize_t abscissasize[1] = {abscissaNslots[i]+1};
        err = H5LTmake_dataset_double(histogram , abscissapath, 1, abscissasize,
                                      abscissavec);
        free(abscissavec);
        if(err){
            H5Gclose(histogram);
            return 1;
        }

        /* Vector length */
        sprintf(temppath, "abscissa_nbin_%02d", i+1);
        err = H5LTmake_dataset_int(histogram, temppath, 1, &dimsize,
                                   &(abscissaNslots[i]));
        if(err){
            H5Gclose(histogram);
            return 1;
        }

        /* Name and unit */
        sprintf(temppath, "name_%02d", i);
        H5LTset_attribute_string(histogram, abscissapath, temppath,
                                 abscissaNames[i]);
        sprintf(temppath, "unit_%02d", i);
        H5LTset_attribute_string(histogram, abscissapath, temppath,
                                 abscissaUnits[i]);

    }

    H5Gclose(histogram);

    return 0;
}
