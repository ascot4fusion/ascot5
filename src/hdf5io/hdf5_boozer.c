/**
 * @file hdf5_boozer.c
 * @brief Module for reading Boozer data from HDF5 file
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../boozer.h"
#include "../consts.h"
#include "hdf5_helpers.h"
#include "hdf5_boozer.h"

/**
 * @brief Initialize Boozer offload data from HDF5 file
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
int hdf5_boozer_init_offload(hid_t f, boozer_offload_data* offload_data,
                             real** offload_array, char* qid) {

    #undef BOOZERPATH
    #define BOOZERPATH "/boozer/Boozer_XXXXXXXXXX/"

    /* Read data the QID corresponds to */
    char path[256];
    hdf5_gen_path(BOOZERPATH, qid, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    /* Read parameters. */
    if( hdf5_read_int(   BOOZERPATH "npsi",       &(offload_data->npsi),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "psimin",    &(offload_data->psi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "psimax",    &(offload_data->psi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(   BOOZERPATH "ntheta", &(offload_data->ntheta),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   BOOZERPATH "nthetag", &(offload_data->nthetag),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(   BOOZERPATH "nrzs", &(offload_data->nrzs),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Size of 1D and 2D input data arrays */
    int nusize      = offload_data->npsi * offload_data->ntheta;
    int thetasize   = offload_data->npsi * offload_data->nthetag;
    int contoursize = offload_data->nrzs;

    /* Allocate offload array */
    *offload_array = (real*) malloc( (nusize + thetasize
                                      + 2 * contoursize) * sizeof(real) );

    /* Read data to offload array */
    if( hdf5_read_double(BOOZERPATH "nu_psitheta",
                         &(*offload_array)[0],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "theta_psithetageom",
                         &(*offload_array)[nusize],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "rs",
                         &(*offload_array)[nusize + thetasize],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "zs",
                         &(*offload_array)[nusize + thetasize
                                           + contoursize],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Initialize the data */
    if( boozer_init_offload(offload_data, offload_array) ) {
        return 1;
    }

    return 0;
}
