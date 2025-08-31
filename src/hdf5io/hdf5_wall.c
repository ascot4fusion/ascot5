/**
 * @file hdf5_wall.c
 * @brief Module for reading wall input from HDF5 file
 *
 * Wall data must be read by calling hdf5_wall_init_offload() contained
 * in this module. This module contains reading routines for all wall data
 * types.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../wall.h"
#include "../wall/wall_2d.h"
#include "../wall/wall_3d.h"
#include "hdf5_wall.h"
#include "hdf5_helpers.h"

#define WPATH /**< Macro that is used to store paths to data groups */

int hdf5_wall_read_2D(hid_t f, wall_2d_data* data, char* qid);
int hdf5_wall_read_3D(hid_t f, wall_3d_data* data, char* qid);

/**
 * @brief Read wall data from HDF5 file
 *
 * This function reads wall data with given qid while also initializing offload
 * data and allocating and filling offload array. The file is opened and closed
 * outside this function.
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to offload data
 * @param offload_array pointer to offload array
 * @param int_offload_array pointer to integer offload array
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_wall_init(hid_t f, wall_data* data, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/wall/wall_2D_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = wall_type_2D;
        err = hdf5_wall_read_2D(f, &(data->w2d), qid);
    }
    hdf5_gen_path("/wall/wall_3D_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = wall_type_3D;
        err = hdf5_wall_read_3D(f, &(data->w3d), qid);
    }
    return err;
}

/**
 * @brief Read 2D wall data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_wall_read_2D(hid_t f, wall_2d_data* data, char* qid) {
    #undef WPATH
    #define WPATH "/wall/wall_2D_XXXXXXXXXX/"

    int nelements;
    if( hdf5_read_int(WPATH "nelements", &nelements,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    real* r = (real*) malloc(nelements * sizeof(real));
    real* z = (real*) malloc(nelements * sizeof(real));
    if( hdf5_read_double(WPATH "r", r, f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(WPATH "z", z, f, qid, __FILE__, __LINE__) ) {return 1;}

    int* flag = (int*) malloc(nelements * sizeof(int));
    if( hdf5_read_int(WPATH "flag", flag,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = wall_2d_init(data, nelements, r, z, flag);
    free(r);
    free(z);
    free(flag);
    return err;
}

/**
 * @brief Read 3D wall data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_wall_read_3D(hid_t f, wall_3d_data* data, char* qid) {
    #undef WPATH
    #define WPATH "/wall/wall_3D_XXXXXXXXXX/"

    int nelements;
    if( hdf5_read_int(WPATH "nelements", &nelements,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate temporary arrays for x1x2x3, y1y2y3, z1z2z3 for each triangle */
    real* x1x2x3 = (real*)malloc(3 * nelements * sizeof(real));
    real* y1y2y3 = (real*)malloc(3 * nelements * sizeof(real));
    real* z1z2z3 = (real*)malloc(3 * nelements * sizeof(real));

    if( hdf5_read_double(WPATH "x1x2x3", x1x2x3,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(WPATH "y1y2y3", y1y2y3,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(WPATH "z1z2z3", z1z2z3,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int* flag = (int*) malloc(nelements * sizeof(int));
    if( hdf5_read_int(WPATH "flag", flag,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    real* vertices;

    int err = wall_3d_init(data, nelements, vertices, flag);
    free(x1x2x3);
    free(y1y2y3);
    free(z1z2z3);
    free(flag);
    return err;
}
