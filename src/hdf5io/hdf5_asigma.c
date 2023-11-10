/**
 * @file hdf5_asigma.c
 * @brief Module for reading atomic reaction cross-section data from HDF5 file
 *
 * Atomic data must be read by calling hdf5_asigma_init_offload()
 * contained in this module. This module contains routines to read atomic
 * reactions data from an ASCOT5 HDF5 file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../asigma.h"
#include "../asigma/asigma_loc.h"
#include "hdf5_helpers.h"
#include "hdf5_asigma.h"

int hdf5_asigma_read_loc(hid_t f, asigma_loc_offload_data* offload_data,
                         real** offload_array, char* qid);

/**
 * @brief Read atomic data from HDF5 file
 *
 * This function reads atomic cross-section (sigma) data with given qid
 * while also initializing offload data and allocating and filling
 * offload array. The file is opened and closed outside this function.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_asigma_init_offload(hid_t f, asigma_offload_data* offload_data,
                             real** offload_array, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/asigma/asigma_loc_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        offload_data->type = asigma_type_loc;
        err = hdf5_asigma_read_loc(f, &(offload_data->asigma_loc),
                                   offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = asigma_init_offload(offload_data, offload_array);
    }

    return err;
}

/**
 * @brief Read atomic sigma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_asigma_read_loc(hid_t f, asigma_loc_offload_data* offload_data,
                         real** offload_array, char* qid) {
    /// @cond
    #undef ASGMPATH
    #define ASGMPATH "/asigma/asigma_loc_XXXXXXXXXX/"
    /// @endcond

    /* Read number of reactions */
    if (hdf5_read_int(ASGMPATH "nreac", &offload_data->N_reac,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read dimensionalities to allow memory alloction for offload array */
    if (hdf5_read_int(ASGMPATH "nenergy", offload_data->N_E,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "ndensity", offload_data->N_n,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "ntemperature", offload_data->N_T,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Helper variables for number of reactions and abscissa dimensions */
    int N_reac = offload_data->N_reac;

    /* Allocate data for the abscissa limits [Emin,Emax,nmin,nmax,Tmin,Tmax]
     * and then for the actual data */
    offload_data->offload_array_length = 6*N_reac;
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        offload_data->offload_array_length +=
              offload_data->N_E[i_reac] * offload_data->N_n[i_reac]
            * offload_data->N_T[i_reac];
    }
    *offload_array = (real*) malloc(offload_data->offload_array_length
                                    *sizeof(real));

    real* E_min     = &(*offload_array)[0*N_reac];
    real* E_max     = &(*offload_array)[1*N_reac];
    real* n_min     = &(*offload_array)[2*N_reac];
    real* n_max     = &(*offload_array)[3*N_reac];
    real* T_min     = &(*offload_array)[4*N_reac];
    real* T_max     = &(*offload_array)[5*N_reac];
    real* sigma     = &(*offload_array)[6*N_reac];

    /* Read reaction identifiers, rest of abscissa parameters and
       reaction data into the offload array*/
    if (hdf5_read_int(ASGMPATH "z1", offload_data->z_1,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "a1", offload_data->a_1,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "z2", offload_data->z_2,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "a2", offload_data->a_2,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "reactype", offload_data->reac_type,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    if (hdf5_read_double(ASGMPATH "energymin", E_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "energymax", E_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "densitymin", n_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "densitymax", n_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "temperaturemin", T_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "temperaturemax", T_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(ASGMPATH "sigma", sigma,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}
