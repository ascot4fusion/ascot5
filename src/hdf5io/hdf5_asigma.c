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

int hdf5_asigma_read_loc(hid_t f, asigma_loc_data* data, char* qid);

/**
 * @brief Read atomic data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to the data
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_asigma_init(hid_t f, asigma_data* data, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */
    hdf5_gen_path("/asigma/asigma_loc_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        data->type = asigma_type_loc;
        err = hdf5_asigma_read_loc(f, &data->asigma_loc, qid);
    }

    return err;
}

/**
 * @brief Read atomic sigma data from HDF5 file
 *
 * @param f HDF5 file from which data is read
 * @param data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_asigma_read_loc(hid_t f, asigma_loc_data* data, char* qid) {
    /// @cond
    #undef ASGMPATH
    #define ASGMPATH "/asigma/asigma_loc_XXXXXXXXXX/"
    /// @endcond

    int nreac;
    if (hdf5_read_int(ASGMPATH "nreac", &nreac,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    int* ne = (int*) malloc( nreac * sizeof(int) );
    int* nn = (int*) malloc( nreac * sizeof(int) );
    int* nT = (int*) malloc( nreac * sizeof(int) );
    int* z1 = (int*) malloc( nreac * sizeof(int) );
    int* a1 = (int*) malloc( nreac * sizeof(int) );
    int* z2 = (int*) malloc( nreac * sizeof(int) );
    int* a2 = (int*) malloc( nreac * sizeof(int) );
    int* reactype = (int*) malloc( nreac * sizeof(int) );
    real* emin = (real*) malloc( nreac * sizeof(real) );
    real* emax = (real*) malloc( nreac * sizeof(real) );
    real* nmin = (real*) malloc( nreac * sizeof(real) );
    real* nmax = (real*) malloc( nreac * sizeof(real) );
    real* Tmin = (real*) malloc( nreac * sizeof(real) );
    real* Tmax = (real*) malloc( nreac * sizeof(real) );
    if (hdf5_read_int(ASGMPATH "nenergy", ne,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "ndensity", nn,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "ntemperature", nT,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "z1", z1,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "a1", a1,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "z2", z2,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "a2", a2,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_int(ASGMPATH "reactype", reactype,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "energymin", emin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "energymax", emax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "densitymin", nmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "densitymax", nmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "temperaturemin", Tmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if (hdf5_read_double(ASGMPATH "temperaturemax", Tmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int nsigmadata = 0;
    for(int i = 0; i < nreac; i++) {
        nsigmadata += ne[i] * nn[i] * nT[i];
    }
    real* sigma = (real*) malloc(nsigmadata * sizeof(real));
    if( hdf5_read_double(ASGMPATH "sigma", sigma,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = asigma_loc_init(data, nreac, z1, a1, z2, a2, reactype,
                              ne, emin, emax, nn, nmin, nmax, nT, Tmin, Tmax,
                              sigma);
    free(ne);
    free(nn);
    free(nT);
    free(z1);
    free(a1);
    free(z2);
    free(a2);
    free(sigma);
    free(reactype);
    free(emin);
    free(emax);
    free(nmin);
    free(nmax);
    free(Tmin);
    free(Tmax);
    return err;
}
