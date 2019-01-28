/**
 * @file hdf5_mhd.c
 * @brief Module for reading MHD data from HDF5 file
 *
 */
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../mhd.h"
#include "hdf5_helpers.h"
#include "hdf5_mhd.h"

/**
 * @brief Initialize MHD offload data from HDF5 file
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is initialized here
 * @param offload_array pointer to offload array which is allocated and
 *                      initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_mhd_init_offload(hid_t f, mhd_offload_data* offload_data,
                          real** offload_array, char* qid) {
    #undef MHDPATH
    #define MHDPATH "/mhd/MHD-XXXXXXXXXX/"

    /* Read data the QID corresponds to */
    char path[256];
    hdf5_gen_path(MHDPATH, qid, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    /* Read parameters. */
    if( hdf5_read_int(MHDPATH "n_modes", &(offload_data->n_modes),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MHDPATH "npsi", &(offload_data->npsi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MHDPATH "nmode", offload_data->nmode,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MHDPATH "mmode", offload_data->mmode,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload array */
    int npsi    = offload_data->npsi;
    int n_modes = offload_data->n_modes;
    offload_data->offload_array_length = n_modes * (2 + 2 * npsi );
    *offload_array = (real*) malloc( offload_data->offload_array_length
                                     * sizeof(real) );

    /* Read data to offload array */
    if( hdf5_read_double(MHDPATH "amplitude_nm", &(*offload_array)[0],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "omega_nm", &(*offload_array)[n_modes],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "alpha_nm", &(*offload_array)[2*n_modes],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "phi_nm", &(*offload_array)[2*n_modes + npsi],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}
