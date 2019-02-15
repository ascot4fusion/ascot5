/**
 * @file hdf5_plasma.c
 * @brief Module for reading plasma input from HDF5 file
 *
 * Plasma data must be read by calling hdf5_plasma_init_offload() contained
 * in this module. This module contains reading routines for all plasma data
 * types.
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../plasma.h"
#include "../plasma/plasma_1D.h"
#include "../plasma/plasma_1DS.h"
#include "../consts.h"
#include "hdf5_helpers.h"
#include "hdf5_plasma.h"

#define PLSPATH /**< Macro that is used to store paths to data groups */

int hdf5_plasma_read_1D(hid_t f, plasma_1D_offload_data* offload_data,
                        real** offload_array, char* qid);
int hdf5_plasma_read_1DS(hid_t f, plasma_1DS_offload_data* offload_data,
                         real** offload_array, char* qid);
/**
 * @brief Read plasma data from HDF5 file
 *
 * This function reads plasma data with given qid while also initializing
 * offload data and allocating and filling offload array. The file is opened
 * and closed outside this function.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_plasma_init_offload(hid_t f, plasma_offload_data* offload_data,
                             real** offload_array, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/plasma/plasma_1D-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        offload_data->type = plasma_type_1D;
        err = hdf5_plasma_read_1D(f, &(offload_data->plasma_1D),
                                  offload_array, qid);
    }

    hdf5_gen_path("/plasma/plasma_1DS-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        offload_data->type = plasma_type_1DS;
        err = hdf5_plasma_read_1DS(f, &(offload_data->plasma_1DS),
                                   offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = plasma_init_offload(offload_data, offload_array);
    }

    return err;
}

/**
 * @brief Read 1D plasma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_plasma_read_1D(hid_t f, plasma_1D_offload_data* offload_data,
                        real** offload_array, char* qid) {
    #undef PLSPATH
    #define PLSPATH "/plasma/plasma_1D-XXXXXXXXXX/"

    /* Read rhogrid size and number of species */
    int n_rho, n_ions;
    if( hdf5_read_int(PLSPATH "n_ions", &n_ions,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "n_rho",  &n_rho,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    offload_data->n_species = n_ions + 1; /* Include electrons */
    offload_data->n_rho = n_rho;


    /* Read species information */
    int temparr[MAX_SPECIES];

    if( hdf5_read_int(PLSPATH "Z_num", temparr,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    for(int i = 0; i < n_ions; i++) {
        offload_data->charge[i+1] = temparr[i] * CONST_E;
    }
    /* Electron charge -1 */
    offload_data->charge[0] = -1 * CONST_E;

    /* Read Amass and calculate mass */
    if( hdf5_read_int(PLSPATH "A_mass", temparr,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    for(int i = 0; i < n_ions; i++) {
        offload_data->mass[i+1] = temparr[i] * CONST_U;
    }
    /* Electron mass */
    offload_data->mass[0] = CONST_M_E;

    /* Allocate space for rhogrid, density (for each species) and
       temperature (for electrons and ions - all ions have same temperature) */
    offload_data->offload_array_length =
        n_rho + 2 * offload_data->n_species*n_rho;
    *offload_array = (real*) malloc(sizeof(real)
                                    * offload_data->offload_array_length);

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* temp_e = &(*offload_array)[n_rho];
    real* temp_i = &(*offload_array)[n_rho*2];
    real* dens_e = &(*offload_array)[n_rho*2
                                     + n_rho*n_ions];
    real* dens_i = &(*offload_array)[n_rho*2
                                     + n_rho*n_ions
                                     + n_rho];

    /* Read rhogrid, densities, and temperatures into allocated array */
    if( hdf5_read_double(PLSPATH "rho", rho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "temp_e", temp_e,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "dens_e", dens_e,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "temp_i", temp_i,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "dens_i", dens_i,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < n_rho; i++) {
        temp_e[i] = temp_e[i] * CONST_E;
        temp_i[i] = temp_i[i] * CONST_E;
    }

    return 0;
}

/**
 * @brief Load plasma data from HDF5 file and prepare parameters
 *
 * This function reads the 1D plasma data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 * @param qid QID of the data
 *
 * @return zero on success
 */
int hdf5_plasma_read_1DS(hid_t f, plasma_1DS_offload_data* offload_data,
                         real** offload_array, char* qid) {

    #undef PLSPATH
    #define PLSPATH "/plasma/plasma_1DS-XXXXXXXXXX/"
    herr_t err = 0;
    int i, j;
    err = H5LTfind_dataset(f, "/plasma/P_1D");

    int n_ions;
    err = H5LTget_attribute_int(f, "/plasma/", "n_ions", &n_ions);
    offload_data->n_species = n_ions + 1; /* Include electrons */

    err = H5LTget_attribute_int(f, "/plasma/P_1D", "n_rho", &(offload_data->n_rho));
    int n_rho = offload_data->n_rho;

    /* Allocate space for rho + temperature for each ion species and electrons
     * + density for each ion species and electrons */
    offload_data->offload_array_length = n_rho + 2*offload_data->n_species*n_rho;
    *offload_array = (real*) malloc(sizeof(real)
                                    * offload_data->offload_array_length);

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* temp_e = &(*offload_array)[n_rho];
    real* temp_i = &(*offload_array)[n_rho*2];
    real* dens_e = &(*offload_array)[n_rho*2
                                     + n_rho*n_ions];
    real* dens_i = &(*offload_array)[n_rho*2
                                     + n_rho*n_ions
                                     + n_rho];

    /* Read Znum and calculate charge */
    int Znum[n_ions];
    err = H5LTread_dataset_int(f,"/plasma/Z_num",Znum);
    /* Electron charge -1 */
    offload_data->charge[0] = -1 * CONST_E;
    for(i = 0; i < n_ions; i++) {
        offload_data->charge[i+1] = Znum[i] * CONST_E;
    }

    /* Read Amass and calculate mass */
    int Amass[n_ions];
    err = H5LTread_dataset_int(f,"/plasma/A_mass", Amass);
    offload_data->mass[0] = CONST_M_E;
    for(i = 0; i < n_ions; i++) {
        offload_data->mass[i+1] = Amass[i] * CONST_U;
    }

    /* Read actual data into array */
    err = H5LTread_dataset_double(f,"/plasma/P_1D/rho", rho);
    err = H5LTread_dataset_double(f,"/plasma/P_1D/temp_e", temp_e);
    for(i = 0; i < n_rho; i++) {
        temp_e[i] = temp_e[i] * CONST_E;
    }
    err = H5LTread_dataset_double(f,"/plasma/P_1D/dens_e", dens_e);

    /* All ions have same temperature */
    err = H5LTread_dataset_double(f,"/plasma/P_1D/temp_i", temp_i);
    for(i = 0; i < n_rho; i++) {
        temp_i[i] = temp_i[i] * CONST_E;
    }

    real temp_dens_i[n_ions*n_rho];
    err = H5LTread_dataset_double(f,"/plasma/P_1D/dens_i", temp_dens_i);
    for(j = 0; j < n_ions; j++) {
        for(i = 0; i < n_rho; i++) {
            dens_i[j*n_rho + i] = temp_dens_i[j*n_rho + i];
        }
    }
    return err;
}
