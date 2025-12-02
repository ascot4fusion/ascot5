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
#include "../plasma/plasma_2D.h"
#include "../plasma/plasma_1Dt.h"
#include "../plasma/plasma_1DS.h"
#include "../consts.h"
#include "hdf5_helpers.h"
#include "hdf5_plasma.h"

#define PLSPATH /**< Macro that is used to store paths to data groups */

int hdf5_plasma_read_1D(hid_t f, plasma_1D_data* data, char* qid);
int hdf5_plasma_read_2D(hid_t f, plasma_2D_data* data, char* qid);
int hdf5_plasma_read_1Dt(hid_t f, plasma_1Dt_data* data, char* qid);
int hdf5_plasma_read_1DS(hid_t f, plasma_1DS_data* data, char* qid);
/**
 * @brief Read plasma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_plasma_init(hid_t f, plasma_data* data, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/plasma/plasma_1D_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = plasma_type_1D;
        err = hdf5_plasma_read_1D(f, &data->plasma_1D, qid);
    }
    hdf5_gen_path("/plasma/plasma_2D_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = plasma_type_2D;
        err = hdf5_plasma_read_2D(f, &data->plasma_2D, qid);
    }
    hdf5_gen_path("/plasma/plasma_1Dt_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = plasma_type_1Dt;
        err = hdf5_plasma_read_1Dt(f, &data->plasma_1Dt, qid);
    }
    hdf5_gen_path("/plasma/plasma_1DS_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = plasma_type_1DS;
        err = hdf5_plasma_read_1DS(f, &data->plasma_1DS, qid);
    }
    return err;
}

/**
 * @brief Read 1D plasma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_plasma_read_1D(hid_t f, plasma_1D_data* data, char* qid) {
    #undef PLSPATH
    #define PLSPATH "/plasma/plasma_1D_XXXXXXXXXX/"

    int n_rho, n_ions;
    if( hdf5_read_int(PLSPATH "nion", &n_ions,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "nrho",  &n_rho,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    int* q = (int*) malloc( n_ions * sizeof(int) );
    int* znum = (int*) malloc( n_ions * sizeof(int) );
    int* anum = (int*) malloc( n_ions * sizeof(int) );
    real* mass = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    real* charge = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    if( hdf5_read_int(PLSPATH "znum", znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "charge", q,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "mass", &mass[1],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    charge[0] = -1 * CONST_E;
    for(int i = 0; i < n_ions; i++) {
        charge[i+1] = q[i] * CONST_E;
    }
    mass[0] = CONST_M_E;
    for(int i = 0; i < n_ions; i++) {
        mass[i+1] *= CONST_U;
    }
    free(q);

    real* Te = (real*) malloc( n_rho*sizeof(real) );
    real* Ti = (real*) malloc( n_rho*sizeof(real) );
    real* ne = (real*) malloc( n_rho*sizeof(real) );
    real* ni = (real*) malloc( n_rho*n_ions*sizeof(real) );
    real* rho = (real*) malloc( n_rho*sizeof(real) );
    real* vtor = (real*) malloc( n_rho*sizeof(real) );
    if( hdf5_read_double(PLSPATH "rho", rho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "vtor", vtor,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "etemperature", Te,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "edensity", ne,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "itemperature", Ti,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "idensity", ni,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < n_rho; i++) {
        Te[i] = Te[i] * CONST_E;
        Ti[i] = Ti[i] * CONST_E;
    }

    int err = plasma_1D_init(data, n_rho, n_ions, rho, anum, znum, mass, charge,
                             Te, Ti, ne, ni, vtor);
    free(Te);
    free(Ti);
    free(ne);
    free(ni);
    free(rho);
    free(znum);
    free(anum);
    free(mass);
    free(charge);
    return err;
}

/**
 * @brief Read 2D plasma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_plasma_read_2D(hid_t f, plasma_2D_data* data, char* qid) {
    #undef PLSPATH
    #define PLSPATH "/plasma/plasma_2D_XXXXXXXXXX/"

    int nr, nz, n_ions;
    real rmin, rmax, zmin, zmax;
    if( hdf5_read_int(PLSPATH "nion", &n_ions,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "nr",  &nr,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "nz",  &nz,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "rmin", &rmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "rmax", &rmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "zmin", &zmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "zmax", &zmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int* q = (int*) malloc( n_ions * sizeof(int) );
    int* znum = (int*) malloc( n_ions * sizeof(int) );
    int* anum = (int*) malloc( n_ions * sizeof(int) );
    real* mass = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    real* charge = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    if( hdf5_read_int(PLSPATH "znum", znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "charge", q,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "mass", &mass[1],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    charge[0] = -1 * CONST_E;
    for(int i = 0; i < n_ions; i++) {
        charge[i+1] = q[i] * CONST_E;
    }
    mass[0] = CONST_M_E;
    for(int i = 0; i < n_ions; i++) {
        mass[i+1] *= CONST_U;
    }
    free(q);

    real* Te = (real*) malloc( nr*nz*sizeof(real) );
    real* Ti = (real*) malloc( nr*nz*sizeof(real) );
    real* ne = (real*) malloc( nr*nz*sizeof(real) );
    real* ni = (real*) malloc( nr*nz*n_ions*sizeof(real) );
    real* vtor = (real*) malloc( nr*nz*sizeof(real) );
    if( hdf5_read_double(PLSPATH "vtor", vtor,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "etemperature", Te,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "edensity", ne,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "itemperature", Ti,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "idensity", ni,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < nr * nz; i++) {
        Te[i] = Te[i] * CONST_E;
        Ti[i] = Ti[i] * CONST_E;
    }

    int err = plasma_2D_init(
        data, nr, nz, n_ions, rmin, rmax, zmin, zmax, anum, znum, mass, charge,
        Te, Ti, ne, ni, vtor);
    free(Te);
    free(Ti);
    free(ne);
    free(ni);
    free(znum);
    free(anum);
    free(mass);
    free(charge);
    return err;
}

/**
 * @brief Read 1Dt plasma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_plasma_read_1Dt(hid_t f, plasma_1Dt_data* data, char* qid) {
    #undef PLSPATH
    #define PLSPATH "/plasma/plasma_1Dt_XXXXXXXXXX/"

    /* Read rhogrid size and number of species */
    int n_rho, n_time, n_ions;
    if( hdf5_read_int(PLSPATH "nion", &n_ions,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "nrho",  &n_rho,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "ntime",  &n_time,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    int* q = (int*) malloc( n_ions * sizeof(int) );
    int* znum = (int*) malloc( n_ions * sizeof(int) );
    int* anum = (int*) malloc( n_ions * sizeof(int) );
    real* mass = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    real* charge = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    if( hdf5_read_int(PLSPATH "znum", znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "charge", q,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "mass", &mass[1],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    charge[0] = -1 * CONST_E;
    for(int i = 0; i < n_ions; i++) {
        charge[i+1] = q[i] * CONST_E;
    }
    mass[0] = CONST_M_E;
    for(int i = 0; i < n_ions; i++) {
        mass[i+1] *= CONST_U;
    }
    free(q);

    real* rho = (real*) malloc( n_rho * sizeof(real) );
    real* time = (real*) malloc( n_time * sizeof(real) );
    if( hdf5_read_double(PLSPATH "rho", rho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "time", time,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    real* vtor = (real*) malloc( n_time*n_rho*sizeof(real) );
    real* Te = (real*) malloc( n_time*n_rho*sizeof(real) );
    real* Ti = (real*) malloc( n_time*n_rho*sizeof(real) );
    real* ne = (real*) malloc( n_time*n_rho*sizeof(real) );
    real* ni = (real*) malloc( n_time*n_ions*n_rho*sizeof(real) );
    if( hdf5_read_double(PLSPATH "etemperature", Te,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "itemperature", Ti,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "vtor", vtor,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i_time = 0; i_time < n_time; i_time++) {
        for(int i_rho = 0; i_rho < n_rho; i_rho++) {
            Te[i_time*n_rho+i_rho] *= CONST_E;
            Ti[i_time*n_rho+i_rho] *= CONST_E;
        }
    }

    if( hdf5_read_double(PLSPATH "edensity", ne,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "idensity", ni,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = plasma_1Dt_init(data, n_rho, n_time, n_ions, rho, time, anum,
                              znum, mass, charge, Te, Ti, ne, ni, vtor);
    free(Te);
    free(Ti);
    free(ne);
    free(ni);
    free(vtor);
    free(znum);
    free(anum);
    free(mass);
    free(charge);
    return err;
}

/**
 * @brief Load plasma data from HDF5 file and prepare parameters
 *
 * @param f hdf5 file identifier
 * @param data pointer to the data struct
 * @param qid QID of the data
 *
 * @return zero on success
 */
int hdf5_plasma_read_1DS(hid_t f, plasma_1DS_data* data, char* qid) {

    #undef PLSPATH
    #define PLSPATH "/plasma/plasma_1DS_XXXXXXXXXX/"

    /* Read rhogrid and number of species */
    int n_rho, n_ions;
    real rhomin, rhomax;
    if( hdf5_read_int(PLSPATH "nion", &n_ions,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "nrho",  &n_rho,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "rhomin",  &rhomin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "rhomax",  &rhomax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int* q = (int*) malloc( n_ions * sizeof(int) );
    int* znum = (int*) malloc( n_ions * sizeof(int) );
    int* anum = (int*) malloc( n_ions * sizeof(int) );
    real* mass = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    real* charge = (real*) malloc( ( n_ions + 1 ) * sizeof(real) );
    if( hdf5_read_int(PLSPATH "znum", znum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(PLSPATH "charge", q,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "mass", &mass[1],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    charge[0] = -1 * CONST_E;
    for(int i = 0; i < n_ions; i++) {
        charge[i+1] = q[i] * CONST_E;
    }
    mass[0] = CONST_M_E;
    for(int i = 0; i < n_ions; i++) {
        mass[i+1] *= CONST_U;
    }
    free(q);

    real* vtor = (real*) malloc( n_rho*sizeof(real) );
    real* Te = (real*) malloc( n_rho*sizeof(real) );
    real* Ti = (real*) malloc( n_rho*sizeof(real) );
    real* ne = (real*) malloc( n_rho*sizeof(real) );
    real* ni = (real*) malloc( n_rho*n_ions*sizeof(real) );
    if( hdf5_read_double(PLSPATH "vtor", vtor,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "etemperature", Te,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "edensity", ne,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "itemperature", Ti,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(PLSPATH "idensity", ni,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i = 0; i < n_rho; i++) {
        Te[i] = Te[i] * CONST_E;
        Ti[i] = Ti[i] * CONST_E;
    }

    int err = plasma_1DS_init(data, n_rho, rhomin, rhomax, n_ions, anum, znum,
                              mass, charge, Te, Ti, ne, ni, vtor);

    free(Te);
    free(Ti);
    free(ne);
    free(ni);
    free(vtor);
    free(znum);
    free(anum);
    free(mass);
    free(charge);
    return err;
}
