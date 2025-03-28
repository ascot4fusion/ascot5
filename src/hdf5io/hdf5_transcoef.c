/**
 * @file hdf5_transcoef.c
 * @brief Module for writing transcport coefficients to a HDF5 file
 */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_helpers.h"
#include "../ascot5.h"
#include "../diag/diag_transcoef.h"
#include "hdf5_transcoef.h"

/**
 * @brief Write transport coefficients to a HDF5 file
 *
 * @param f hdf5 file
 * @param path path to group which is created here and where the data is stored
 * @param data transport coefficient diagnostics offload data
 * @param coefarr array storing the coefficient data [id, K, D]
 *
 * @return zero on success
 */
int hdf5_transcoef_write(hid_t f, char* path, diag_transcoef_data* data) {
    hid_t group = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(group < 0) {
        return 1;
    }

    /* Find number of (valid) data points in the coefarray */
    int arrlen = data->Nmrk;
    int datasize = 0;
    integer* mask = malloc(arrlen*sizeof(real));
    for(integer i=0; i < arrlen; i++) {
        mask[i] = data->id[i] > 0;
        if(mask[i]) {
            datasize++;
        }
    }

    /* Write IDs*/
    integer j = 0;
    integer* idarr = (integer*) malloc(datasize * sizeof(integer));
    for(integer i = 0; i < arrlen; i++) {
        if(mask[i]) {
            idarr[j] = data->id[i];
            j++;
        }
    }
    hdf5_write_extendible_dataset_long(group, "ids", datasize, idarr);
    free(idarr);

    /* Write K and D */
    j = 0;
    real* dataarr = (real*) malloc(datasize * sizeof(real));
    for(integer i = 0; i < arrlen; i++) {
        if(mask[i]) {
            dataarr[j] = data->Kcoef[i];
            j++;
        }
    }
    hdf5_write_extendible_dataset_double(group, "k", datasize, dataarr);

    j = 0;
    for(integer i = 0; i < arrlen; i++) {
        if(mask[i]) {
            dataarr[j] = data->Dcoef[i];
            j++;
        }
    }
    hdf5_write_extendible_dataset_double(group, "d", datasize, dataarr);
    free(dataarr);
    free(mask);

    /* Write units */
    H5LTset_attribute_string(group, "ids", "unit", "1");
    if( data->recordrho ) {
        H5LTset_attribute_string(group, "k", "unit", "1/s");
        H5LTset_attribute_string(group, "d", "unit", "1/s");
    }
    else {
        H5LTset_attribute_string(group, "k", "unit", "m/s");
        H5LTset_attribute_string(group, "d", "unit", "m^2/s");
    }

    H5Gclose (group);

    return 0;
}
