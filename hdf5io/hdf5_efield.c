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
#include "../Efield/E_1D.h"
#include "../Efield/E_1DS.h"
#include "../print.h"
#include "hdf5_helpers.h"
#include "hdf5_efield.h"

int hdf5_efield_read_1D(hid_t f, E_1D_offload_data* offload_data,
                        real** offload_array, char* qid);
int hdf5_efield_read_1DS(hid_t f, E_1DS_offload_data* offload_data,
                         real** offload_array, char* qid);
int hdf5_efield_read_TC(hid_t f, E_TC_offload_data* offload_data,
                        real** offload_array, char* qid);
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
 * @param offload_data pointer to offload data struct which is initialized here
 * @param offload_array pointer to offload array which is allocated and
 *                      initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data,
                             real** offload_array, char* qid) {
    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/efield/E_TC-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = E_field_type_TC;
        err = hdf5_efield_read_TC(f, &(offload_data->ETC),
                                  offload_array, qid);
    }

    hdf5_gen_path("/efield/E_1DS-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        offload_data->type = E_field_type_1DS;
        err = hdf5_efield_read_1DS(f, &(offload_data->E1DS),
                                   offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = E_field_init_offload(offload_data, offload_array);
    }

    return err;
}

/**
 * @brief Read E1D electric field data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_efield_read_1D(hid_t f, E_1D_offload_data* offload_data,
                        real** offload_array, char* qid) {

    herr_t err;
    char path[256];

    err = H5LTget_attribute_int(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "n_rho", &(offload_data->n_rho));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "rho_min", &(offload_data->rho_min));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "rho_max", &(offload_data->rho_max));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    
    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/rho", qid, path), rho);
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/dV_drho", qid, path), dV);
    if(err) {print_err("Error: Failed to read dataset."); return 1;}

    /* Effective minor radius */
    real r_eff;
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "r_eff", &(r_eff));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    
    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }

    return err;
}

/**
 * @brief Read E1DS electric field data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_efield_read_1DS(hid_t f, E_1DS_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_1DS-XXXXXXXXXX/"

    if( hdf5_read_int(EPATH "n_rho", &(offload_data->n_rho),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "rho_min", &(offload_data->rho_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "rho_max", &(offload_data->rho_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];

    if( hdf5_read_double(EPATH "rho", rho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "dV_drho", dV,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Effective minor radius */
    real r_eff;
    if( hdf5_read_double(EPATH "r_eff", &(r_eff),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }

    return 0;
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
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array but no data is stored there so
 *                      it is not allocated and NULL pointer is returned instead
 * @param qid QID of the B_TC field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_efield_read_TC(hid_t f, E_TC_offload_data* offload_data,
                        real** offload_array, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_TC-XXXXXXXXXX/"

    *offload_array = NULL;

    if( hdf5_read_double(EPATH "Exyz", offload_data->Exyz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}
