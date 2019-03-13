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
    #define MHDPATH "/mhd/MHD_XXXXXXXXXX/"

    /* Read data the QID corresponds to */
    char path[256];
    hdf5_gen_path(MHDPATH, qid, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    /* Read parameters. */
    if( hdf5_read_int(   MHDPATH "nmode",      &(offload_data->n_modes),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "npsi",       &(offload_data->npsi),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "psimin",     &(offload_data->psi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "psimax",     &(offload_data->psi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "ntime",      &(offload_data->ntime),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "tmin",       &(offload_data->t_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "tmax",       &(offload_data->t_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "nmodes",       offload_data->nmode,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "mmodes",       offload_data->mmode,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "amplitude",    offload_data->amplitude_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "omega",        offload_data->omega_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload array */
    int datasize = offload_data->npsi * offload_data->ntime
                   * offload_data->n_modes;
    *offload_array = (real*) malloc( 2 * datasize * sizeof(real) );

    /* Read data to offload array */
    if( hdf5_read_double(MHDPATH "alpha", &(*offload_array)[0],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "phi",   &(*offload_array)[datasize],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Initialize the data */
    if( mhd_init_offload(offload_data, offload_array) ) {
        return 1;
    }

    return 0;
}
