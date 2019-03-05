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
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../neutral.h"
#include "../neutral/N0_3D.h"
#include "../consts.h"
#include "../math.h"
#include "hdf5_neutral.h"
#include "hdf5_helpers.h"

#define NPATH /**< Macro that is used to store paths to data groups */

int hdf5_neutral_read_3D(hid_t f, N0_3D_offload_data* offload_data,
                         real** offload_array, char* qid);

/**
 * @brief Initialize neutral data from HDF5 file
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
int hdf5_neutral_init_offload(hid_t f, neutral_offload_data* offload_data,
                              real** offload_array, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/neutral/N0_3D_XXXXXXXXXX", qid, path);
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
 * @param qid QID of the data that is to be read
 *
 * @return zero on success
 */
int hdf5_neutral_read_3D(hid_t f, N0_3D_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef NPATH
    #define NPATH "/neutral/N0_3D_XXXXXXXXXX/"

    /* Read and initialize Rpz-grid */
    if( hdf5_read_int(NPATH "nr", &(offload_data->n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "nphi", &(offload_data->n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "nz", &(offload_data->n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "rmin", &(offload_data->r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "rmax", &(offload_data->r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phimin", &(offload_data->phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phimax", &(offload_data->phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "zmin", &(offload_data->z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "zmax", &(offload_data->z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Convert to radians */
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);

    /* Read n_species, anum, znum and distribution type */
    if( hdf5_read_int(NPATH "n_species", &(offload_data->n_species),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "anum", offload_data->anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "znum", offload_data->znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "maxwellian", offload_data->maxwellian,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    int N0_size = offload_data->n_r * offload_data->n_phi * offload_data->n_z;
    int T0_size = offload_data->n_r * offload_data->n_phi * offload_data->n_z;

    *offload_array = (real*) malloc(offload_data->n_species
                                    * (N0_size + T0_size)
                                    * sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* n0 = &(*offload_array)[0];
    real* t0 = &(*offload_array)[offload_data->n_species * N0_size];

    /* Read the neutral density and temperature */
    if( hdf5_read_double(NPATH "density", n0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "temperature", t0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < offload_data->n_species * T0_size; i++) {
        t0[i] = t0[i] * CONST_E;
    }

    return 0;
}
