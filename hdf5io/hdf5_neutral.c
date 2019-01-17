/**
 * @file hdf5_neutral.c
 * @brief Module for reading neutral data from HDF5 file
 *
 * Neutral data  must be read by calling hdf5_neutral_init_offload() contained
 * in this module. This module contains reading routines for all neutral data
 * types.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ascot5.h"
#include "../neutral.h"
#include "../neutral/N0_3D.h"
#include "../consts.h"
#include "../math.h"
#include "hdf5_neutral.h"
#include <hdf5.h>
#include "hdf5_helpers.h"
#include <hdf5_hl.h>

int hdf5_neutral_read_3D(hid_t f, N0_3D_offload_data* offload_data,
                         real** offload_array, char* qid);

int hdf5_neutral_read_ST(hid_t f, N0_ST_offload_data* offload_data,
                         real** offload_array, char* qid);


int hdf5_neutral_init_offload(hid_t f, neutral_offload_data* offload_data,
                              real** offload_array, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/neutral/N0_3D-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = neutral_type_3D;
        err = hdf5_neutral_read_3D(f, &(offload_data->N03D),
                                   offload_array, qid);
    }

    hdf5_gen_path("/neutral/N0_3D-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = neutral_type_3D;
        err = hdf5_neutral_read_3D(f, &(offload_data->N03D),
                                   offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = neutral_init_offload(offload_data, offload_array);
    }

    return err;
}

/**
 * @brief Load neutral data from HDF5 file and prepare parameters
 *
 * This function reads the 3D neutral data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 */
int hdf5_neutral_read_3D(hid_t f, N0_3D_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef NPATH
    #define NPATH "/neutral/N0_3D-XXXXXXXXXX/"

    /* Read and initialize Rpz-grid */
    if( hdf5_read_int(NPATH "n_r", &(offload_data->n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "n_phi", &(offload_data->n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "n_z", &(offload_data->n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "r_min", &(offload_data->r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "r_max", &(offload_data->r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phi_min", &(offload_data->phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phi_max", &(offload_data->phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "z_min", &(offload_data->z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "z_max", &(offload_data->z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);

    *offload_array = (real*) malloc(offload_data->n_r * offload_data->n_phi
                                    * offload_data->n_z * sizeof(real));

    /* Read the neutral density */
    if( hdf5_read_double(NPATH "n0", *offload_array,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    return 0;
}

/**
 * @brief Load neutral data from HDF5 file and prepare parameters
 *
 * This function reads the stellarator neutral data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 */
int hdf5_neutral_read_ST(hid_t f, N0_ST_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef NPATH
    #define NPATH "/neutral/N0_ST-XXXXXXXXXX/"

    if( hdf5_read_int(NPATH "toroidalPeriods", &(offload_data->periods),
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read and initialize Rpz-grid */
    if( hdf5_read_int(NPATH "n_r", &(offload_data->n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "n_phi", &(offload_data->n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "n_z", &(offload_data->n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "r_min", &(offload_data->r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "r_max", &(offload_data->r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phi_min", &(offload_data->phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phi_max", &(offload_data->phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "z_min", &(offload_data->z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "z_max", &(offload_data->z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);

    *offload_array = (real*) malloc(offload_data->n_r * offload_data->n_phi
                                    * offload_data->n_z * sizeof(real));

    /* Read the neutral density */
    if( hdf5_read_double(NPATH "n0", *offload_array,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}
