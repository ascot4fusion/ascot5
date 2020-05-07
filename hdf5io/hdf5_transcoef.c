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
 * @brief Write transcport coefficients to a HDF5 file
 *
 * @param f hdf5 file
 * @param qid qid of the results group
 * @param data transport coefficient diagnostics offload data
 * @param coefdata array storing the coefficient data [id, K, D]
 *
 * @return zero on success
 */
int hdf5_transcoef_write(hid_t f, char* qid, diag_transcoef_offload_data* data,
                         real* coefarr) {
    char path[256];
    hdf5_generate_qid_path("/results/run_XXXXXXXXXX/", qid, path);
    strcat(path, "transcoef");

    hid_t group = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if(group < 0) {
        return 1;
    }

    /* Find number of (valid) data points in the coefarray */
    int arrlen = data->Nmrk;
    int datasize = 0;
    integer* mask = malloc(arrlen*sizeof(real));
    for(integer i=0; i < arrlen; i++) {
        mask[i] = coefarr[i] > 0;
        if(mask[i]) {
            datasize++;
        }
    }

    /* Write IDs*/
    integer j = 0;
    integer* idarr = (integer*) malloc(datasize * sizeof(integer));
    for(integer i = 0; i < arrlen; i++) {
        if(mask[i]) {
            idarr[j] = coefarr[i];
            j++;
        }
    }
    hdf5_write_extendible_dataset_long(group, "id", datasize, idarr);
    free(idarr);

    /* Write K and D */
    j = 0;
    real* dataarr = (real*) malloc(datasize * sizeof(real));
    for(integer i = 0; i < arrlen; i++) {
        if(mask[i]) {
            dataarr[j] = coefarr[i + arrlen];
            j++;
        }
    }
    hdf5_write_extendible_dataset_double(group, "K", datasize, dataarr);

    j = 0;
    for(integer i = 0; i < arrlen; i++) {
        if(mask[i]) {
            dataarr[j] = coefarr[i + 2*arrlen];
            j++;
        }
    }
    hdf5_write_extendible_dataset_double(group, "D", datasize, dataarr);
    free(dataarr);
    free(mask);

    H5Gclose (group);

    return 0;
}
