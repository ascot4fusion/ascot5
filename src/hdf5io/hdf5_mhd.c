/**
 * @file hdf5_mhd.c
 * @brief Module for reading MHD data from HDF5 file
 *
 */
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../mhd.h"
#include "../mhd/mhd_stat.h"
#include "../mhd/mhd_nonstat.h"
#include "hdf5_helpers.h"
#include "hdf5_mhd.h"

#define MHDPATH /**< Macro that is used to store paths to data groups */

int hdf5_mhd_read_stat(hid_t f, mhd_stat_data* data, char* qid);
int hdf5_mhd_read_nonstat(hid_t f, mhd_nonstat_data* data, char* qid);

/**
 * @brief Initialize MHD data from HDF5 file
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param data pointer to the data struct which is initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_mhd_init(hid_t f, mhd_data* data, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/mhd/MHD_STAT_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = mhd_type_stat;
        err = hdf5_mhd_read_stat(f, &data->stat, qid);
    }
    hdf5_gen_path("/mhd/MHD_NONSTAT_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = mhd_type_nonstat;
        err = hdf5_mhd_read_nonstat(f, &data->nonstat, qid);
    }
    return err;
}

/**
 * @brief Read stationary MHD data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_mhd_read_stat(hid_t f, mhd_stat_data* data, char* qid) {
    #undef MHDPATH
    #define MHDPATH "/mhd/MHD_STAT_XXXXXXXXXX/"

    int nmode, nrho;
    real rhomin, rhomax;
    if( hdf5_read_int(   MHDPATH "nmode", &nmode,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "nrho", &nrho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "rhomin", &rhomin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "rhomax", &rhomax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int* moden = (int*) malloc( nmode * sizeof(int) );
    int* modem = (int*) malloc( nmode * sizeof(int) );
    real* omega_nm = (real*) malloc( nmode * sizeof(real) );
    real* phase_nm = (real*) malloc( nmode * sizeof(real) );
    real* amplitude_nm = (real*)malloc( nmode * sizeof(real) );
    if( hdf5_read_int(   MHDPATH "nmodes", moden,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "mmodes", modem,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "omega", omega_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "phase", phase_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "amplitude", amplitude_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload array */
    real* phi = (real*) malloc( nrho * nmode * sizeof(real) );
    real* alpha = (real*) malloc( nrho * nmode * sizeof(real) );
    if( hdf5_read_double(MHDPATH "alpha", alpha,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "phi", phi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int err = mhd_stat_init(data, nmode, nrho, rhomin, rhomax, moden, modem,
                            amplitude_nm, omega_nm, phase_nm, alpha, phi);
    return err;
}

/**
 * @brief Read nonstationary MHD data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_mhd_read_nonstat(hid_t f, mhd_nonstat_data* data, char* qid) {
    #undef MHDPATH
    #define MHDPATH "/mhd/MHD_NONSTAT_XXXXXXXXXX/"

    int nmode, nrho, ntime;
    real rhomin, rhomax, tmin, tmax;
    if( hdf5_read_int(   MHDPATH "nmode", &nmode,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "nrho", &nrho,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "rhomin", &rhomin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "rhomax", &rhomax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "ntime", &ntime,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "tmin", &tmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "tmax", &tmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int* moden = (int*) malloc( nmode * sizeof(int) );
    int* modem = (int*) malloc( nmode * sizeof(int) );
    real* omega_nm = (real*) malloc( nmode * sizeof(real) );
    real* phase_nm = (real*) malloc( nmode * sizeof(real) );
    real* amplitude_nm = (real*)malloc( nmode * sizeof(real) );
    if( hdf5_read_int(   MHDPATH "nmodes", moden,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(   MHDPATH "mmodes", modem,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "omega", omega_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "phase", phase_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "amplitude", amplitude_nm,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload array */
    real* phi = (real*) malloc( nrho * ntime * nmode * sizeof(real) );
    real* alpha = (real*) malloc( nrho * ntime * nmode * sizeof(real) );
    if( hdf5_read_double(MHDPATH "alpha", alpha,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MHDPATH "phi", phi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int err = mhd_nonstat_init(data, nmode, nrho, ntime, rhomin, rhomax,
                               tmin, tmax, moden, modem, amplitude_nm,
                               omega_nm, phase_nm, alpha, phi);
    return err;
}
