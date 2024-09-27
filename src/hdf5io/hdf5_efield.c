/**
 * @file hdf5_efield.c
 * @brief Module for reading electric field data from HDF5 file
 *
 * Electric field  must be read by calling hdf5_efield_init_offload() contained
 * in this module. This module contains reading routines for all electric field
 * types.
 *
 * Note: the routines within this module that read the data from HDF5 file
 * may use the offload data struct and offload arrays as (temporary) storage.
 * However, the actual initialization is done at the specific
 * E_field_init_offload() function the electric field data corresponds to. Check
 * from that function what the offload data and the offload array are expected
 * to contain. As a rule of thumb, the reading routines here should only read
 * the data and maybe do some trivial computations but nothing complicated.
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../E_field.h"
#include "../Efield/E_TC.h"
#include "../Efield/E_1DS.h"
#include "../print.h"
#include "hdf5_helpers.h"
#include "hdf5_efield.h"

#define EPATH /**< Macro that is used to store paths to data groups */

int hdf5_efield_read_1DS(hid_t f, E_1DS_data* data, char* qid);
int hdf5_efield_read_TC(hid_t f, E_TC_data* data, char* qid);

/**
 * @brief Read electric field data from HDF5 file
 *
 * This function reads electric field data with given qid while also
 * initializing offload data and allocating and filling offload array. The file
 * is opened and closed outside this function.
 *
 * The electric field data is stored under /efield/ group in ASCOT5 HDF5 file.
 * Several electric fields of same or different type maybe stored there in
 * their respective groups as long as the group name contains QID as an
 * identifier.
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param data pointer to the data struct which is initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_efield_init(hid_t f, E_field_data* data, char* qid) {
    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/efield/E_TC_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = E_field_type_TC;
        err = hdf5_efield_read_TC(f, &data->ETC, qid);
    }
    hdf5_gen_path("/efield/E_1DS_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = E_field_type_1DS;
        err = hdf5_efield_read_1DS(f, &data->E1DS, qid);
    }
    return err;
}

/**
 * @brief Read E1DS electric field data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_efield_read_1DS(hid_t f, E_1DS_data* data, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_1DS_XXXXXXXXXX/"

    int nrho;
    real rhomin, rhomax, reff;
    if( hdf5_read_int(EPATH "nrho", &nrho,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "rhomin", &rhomin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "rhomax", &rhomax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "reff", &reff,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    real* dvdrho = (real*) malloc( nrho*sizeof(real) );
    if( hdf5_read_double(EPATH "dvdrho", dvdrho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = E_1DS_init(data, nrho, rhomin, rhomax, reff, dvdrho);
    free(dvdrho);
    return err;
}

/**
 * @brief Read magnetic field data of type E_TC
 *
 * The E_TC data is stored in HDF5 file under the group
 * /efield/E_TC-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - double Exyz Electric field values [E_x, E_y, E_z] [V/m]
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param data pointer to the data struct which is allocated here
 * @param qid QID of the B_TC field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_efield_read_TC(hid_t f, E_TC_data* data, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_TC_XXXXXXXXXX/"

    real exyz[3];
    if( hdf5_read_double(EPATH "exyz", exyz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = E_TC_init(data, exyz);
    return err;
}
