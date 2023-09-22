/**
 * @file hdf5_nbi.c
 * @brief Module for reading NBI data from HDF5 file
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h> //?
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../nbi.h"
#include "../consts.h"
#include "hdf5_helpers.h"
#include "hdf5_nbi.h"

/**
 * @brief Initialize NBI offload data from HDF5 file
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
int hdf5_nbi_init_offload(hid_t f, nbi_offload_data* offload_data,
                          real** offload_array) {
    herr_t err0;

    // TODO
    #if VERBOSE > 0
        printf("Reading NBI input from the HDF5 file...\n");
    #endif

    err0 = hdf5_find_group(f, "/nbi/");
    if(err0 < 0) {
        return -1;//TODO?
    }

    char active[11];
    err0 = H5LTget_attribute_string(f, "/nbi/", "active", active);
    if(err0 < 0) {
        return -1;//TODO?
    }
    active[10] = '\0';

    // TODO
    #if VERBOSE > 0
        printf("Active qid is %s\n", active);
    #endif

    char path[256]; // Storage array required for hdf5_gen_path() calls

    /* Read data the QID corresponds to */

    hdf5_generate_qid_path("/nbi/nbi_XXXXXXXXXX/", active, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    #undef NBIPATH
    #define NBIPATH "/nbi/nbi_XXXXXXXXXX/"

    if( hdf5_read_int(NBIPATH "ninj", &(offload_data->ninj),
                      f, active, __FILE__, __LINE__) ) {return 1;}

    /* Loop over injectors twice: first time we initialize offload data */
    offload_data->offload_array_length = 0;
    for(int i = 0; i < offload_data->ninj; i++) {
        sprintf(path, NBIPATH "inj%d/%s", i+1, "ids");
        if( hdf5_read_int(path, &(offload_data->id[i]),
                          f, active, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "nbeamlet");
        if( hdf5_read_int(path, &(offload_data->n_beamlet[i]),
                          f, active, __FILE__, __LINE__) ) {return 1;}

        offload_data->offload_array_length += 6*offload_data->n_beamlet[i];

        sprintf(path, NBIPATH "inj%d/%s", i+1, "power");
        if( hdf5_read_double(path, &(offload_data->power[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "energy");
        if( hdf5_read_double(path, &(offload_data->energy[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "efrac");
        if( hdf5_read_double(path, &(offload_data->efrac[i*3]),
                             f, active, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "divh");
        if( hdf5_read_double(path, &(offload_data->div_h[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divv");
        if( hdf5_read_double(path, &(offload_data->div_v[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divhalofrac");
        if( hdf5_read_double(path, &(offload_data->div_halo_frac[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divhaloh");
        if( hdf5_read_double(path, &(offload_data->div_halo_h[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divhalov");
        if( hdf5_read_double(path, &(offload_data->div_halo_v[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "anum");
        if( hdf5_read_int(path, &(offload_data->anum[i]),
                          f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "znum");
        if( hdf5_read_int(path, &(offload_data->znum[i]),
                          f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "mass");
        if( hdf5_read_double(path, &(offload_data->mass[i]),
                             f, active, __FILE__, __LINE__) ) {return 1;}

        /* Conver to SI */
        offload_data->energy[i] *= CONST_E;
        offload_data->mass[i]   *= CONST_U;
    }

    /* In second loop we now initialize the offload array which we are
     * able to allocate now */
    *offload_array = (real*)malloc(
        offload_data->offload_array_length * sizeof(real) );
    int idx = 0;
    for(int i = 0; i < offload_data->ninj; i++) {
        int nbeamlet = offload_data->n_beamlet[i];
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletx");
        if( hdf5_read_double(path, &(*offload_array)[idx + 0*nbeamlet],
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamlety");
        if( hdf5_read_double(path, &(*offload_array)[idx + 1*nbeamlet],
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletz");
        if( hdf5_read_double(path, &(*offload_array)[idx + 2*nbeamlet],
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletdx");
        if( hdf5_read_double(path, &(*offload_array)[idx + 3*nbeamlet],
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletdy");
        if( hdf5_read_double(path, &(*offload_array)[idx + 4*nbeamlet],
                             f, active, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletdz");
        if( hdf5_read_double(path, &(*offload_array)[idx + 5*nbeamlet],
                             f, active, __FILE__, __LINE__) ) {return 1;}
        idx += 6*nbeamlet;
    }

    /* Initialize the data */
    if( nbi_init_offload(offload_data, offload_array) ) {
        return 1;
    }
    return 0;
}
