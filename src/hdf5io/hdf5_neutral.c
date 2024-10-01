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
#include "../neutral/N0_1D.h"
#include "../neutral/N0_3D.h"
#include "../consts.h"
#include "../math.h"
#include "hdf5_neutral.h"
#include "hdf5_helpers.h"

#define NPATH /**< Macro that is used to store paths to data groups */

int hdf5_neutral_init_1D(hid_t f, N0_1D_data* data, char* qid);
int hdf5_neutral_init_3D(hid_t f, N0_3D_data* data, char* qid);

/**
 * @brief Initialize neutral data from HDF5 file
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param data pointer to the data struct which is initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_neutral_init(hid_t f, neutral_data* data, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/neutral/N0_1D_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = neutral_type_1D;
        err = hdf5_neutral_init_1D(f, &(data->N01D), qid);
    }
    hdf5_gen_path("/neutral/N0_3D_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = neutral_type_3D;
        err = hdf5_neutral_init_3D(f, &(data->N03D), qid);
    }
    return err;
}

/**
 * @brief Load 1D neutral data from HDF5 file and initialize it
 *
 * @param f hdf5 file identifier
 * @param data pointer to the data struct
 * @param qid QID of the data that is to be read
 *
 * @return zero on success
 */
int hdf5_neutral_init_1D(hid_t f, N0_1D_data* data, char* qid) {
    #undef NPATH
    #define NPATH "/neutral/N0_1D_XXXXXXXXXX/"

    /* Read and initialize rho coordinate */
    int n_rho;
    real rho_min, rho_max;
    if( hdf5_read_int(NPATH "nrho", &n_rho,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "rhomin", &rho_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "rhomax", &rho_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read n_species, anum, znum and distribution type */
    int n_species;
    if( hdf5_read_int(NPATH "nspecies", &n_species,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    int* anum = (int*) malloc(n_species * sizeof(int));
    int* znum = (int*) malloc(n_species * sizeof(int));
    int* maxwellian = (int*) malloc(n_species * sizeof(int));
    if( hdf5_read_int(NPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "znum", znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "maxwellian", maxwellian,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the neutral density and temperature */
    real* density = (real*) malloc(n_species*n_rho*sizeof(real));
    real* temperature = (real*) malloc(n_species*n_rho*sizeof(real));
    if( hdf5_read_double(NPATH "density", density,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "temperature", temperature,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < n_species * n_rho; i++) {
        temperature[i] = temperature[i] * CONST_E;
    }

    int err = N0_1D_init(data, n_rho, rho_min, rho_max, n_species, anum, znum,
                         maxwellian, density, temperature);
    free(anum);
    free(znum);
    free(maxwellian);
    free(density);
    free(temperature);

    return err;
}

/**
 * @brief Load 3D neutral data from HDF5 file and initialize it
 *
 * @param f hdf5 file identifier
 * @param data pointer to the data struct
 * @param qid QID of the data that is to be read
 *
 * @return zero on success
 */
int hdf5_neutral_init_3D(hid_t f, N0_3D_data* data, char* qid) {
    #undef NPATH
    #define NPATH "/neutral/N0_3D_XXXXXXXXXX/"

    /* Read and initialize Rpz-grid */
    int n_r, n_phi, n_z;
    real r_min, r_max, phi_min, phi_max, z_min, z_max;
    if( hdf5_read_int(NPATH "nr", &n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "nphi", &n_phi,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "nz", &n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "rmin", &r_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "rmax", &r_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phimin", &phi_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "phimax", &phi_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "zmin", &z_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "zmax", &z_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Convert to radians */
    phi_max = math_deg2rad(phi_max);
    phi_min = math_deg2rad(phi_min);

    /* Read n_species, anum, znum and distribution type */
    int n_species;
    if( hdf5_read_int(NPATH "nspecies", &n_species,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    int* anum = (int*) malloc(n_species * sizeof(int));
    int* znum = (int*) malloc(n_species * sizeof(int));
    int* maxwellian = (int*) malloc(n_species * sizeof(int));
    if( hdf5_read_int(NPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "znum", znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(NPATH "maxwellian", maxwellian,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the neutral density and temperature */
    real* density = (real*) malloc(n_species*n_r*n_phi*n_z*sizeof(real));
    real* temperature = (real*) malloc(n_species*n_r*n_phi*n_z*sizeof(real));
    if( hdf5_read_double(NPATH "density", density,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(NPATH "temperature", temperature,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < n_species*n_r*n_phi*n_z; i++) {
        temperature[i] = temperature[i] * CONST_E;
    }

    int err = N0_3D_init(data, n_r, r_min, r_max, n_phi, phi_min, phi_max,
                         n_z, z_min, z_max, n_species, anum, znum, maxwellian,
                         density, temperature);
    free(anum);
    free(znum);
    free(maxwellian);
    free(density);
    free(temperature);
    return err;
}
