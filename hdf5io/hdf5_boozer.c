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
    #define BOOZERPATH "/boozer/Boozer-XXXXXXXXXX/"

    /* Read data the QID corresponds to */
    char path[256];
    hdf5_gen_path(BOOZERPATH, qid, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    /* Read parameters. */
    if( hdf5_read_int(BOOZERPATH "nR", &(offload_data->nR),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "R_min", &(offload_data->R_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "R_max", &(offload_data->R_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BOOZERPATH "nz", &(offload_data->nz),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "z_min", &(offload_data->z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "z_max", &(offload_data->z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BOOZERPATH "npsi", &(offload_data->npsi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "psi_min", &(offload_data->psi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "psi_max", &(offload_data->psi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BOOZERPATH "ntheta_geo", &(offload_data->ntheta_geo),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BOOZERPATH "ntheta_bzr", &(offload_data->ntheta_bzr),
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload array */
    int npsi       = offload_data->npsi;
    int ntheta_geo = offload_data->ntheta_geo;
    int ntheta_bzr = offload_data->ntheta_bzr;
    int Rzgridsize = offload_data->nR * offload_data->nz;
    offload_data->offload_array_length =
        npsi * (3 + ntheta_geo + 2*ntheta_bzr) + Rzgridsize;
    *offload_array = (real*) malloc( offload_data->offload_array_length
                                     * sizeof(real) );

    /* Read data to offload array */
    if( hdf5_read_double(BOOZERPATH "g", &(*offload_array)[0],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "q", &(*offload_array)[npsi],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "I", &(*offload_array)[2*npsi],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "delta", &(*offload_array)[3*npsi],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "nu",
                         &(*offload_array)[npsi*(3+ntheta_bzr)],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "theta_bzr",
                         &(*offload_array)[npsi*(3+2*ntheta_bzr)],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "theta_geo",
                         &(*offload_array)[npsi*(3+2*ntheta_bzr+ntheta_geo)],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Initialize the data */
    if( boozer_init_offload(offload_data, offload_array) ) {
        return 1;
    }

    return 0;
}
