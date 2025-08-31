/**
 * @file hdf5_boozer.c
 * @brief Module for reading Boozer data from HDF5 file
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../boozer.h"
#include "../consts.h"
#include "hdf5_helpers.h"
#include "hdf5_boozer.h"

/**
 * @brief Initialize Boozer data from HDF5 file
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param data pointer to the data struct which is initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_boozer_init(hid_t f, boozer_data* data, char* qid) {

    /// @cond
    #undef BOOZERPATH
    #define BOOZERPATH "/boozer/Boozer_XXXXXXXXXX/"
    /// @endcond

    /* Read data the QID corresponds to */
    char path[256];
    hdf5_gen_path(BOOZERPATH, qid, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    /* Read parameters. */
    int ntheta, nthetag, nrzs, npsi;
    real psimin, psimax;
    if( hdf5_read_int(BOOZERPATH "ntheta", &ntheta,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BOOZERPATH "nthetag", &nthetag,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BOOZERPATH "nrzs", &nrzs,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BOOZERPATH "npsi", &npsi,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "psimin", &psimin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "psimax", &psimax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read data to offload array */
    real* rs = (real*) malloc( npsi * nrzs * sizeof(real) );
    real* zs = (real*) malloc( npsi * nrzs * sizeof(real) );
    real* nu = (real*) malloc( npsi * ntheta * sizeof(real) );
    real* theta = (real*) malloc( npsi * nthetag * sizeof(real) );
    if( hdf5_read_double(BOOZERPATH "nu_psitheta", nu,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "theta_psithetageom", theta,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "rs", rs,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BOOZERPATH "zs", zs,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int err = boozer_init(data, npsi, psimin, psimax, ntheta, nthetag, 4,
                          nu, theta, nrzs, rs, zs);
    free(rs);
    free(zs);
    free(nu);
    free(theta);
    return err;
}
