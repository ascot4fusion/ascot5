/**
 * @file hdf5_nbi.c
 * @brief Module for reading NBI data from HDF5 file
 */
#include <stdio.h>
#include <stdlib.h>
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
 * @param data pointer to the data struct which is initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_nbi_init(hid_t f, nbi_data* data, char* qid) {

    char path[256]; // Storage array required for hdf5_gen_path() calls

    /* Read data the QID corresponds to */

    hdf5_generate_qid_path("/nbi/nbi_XXXXXXXXXX/", qid, path);
    if( hdf5_find_group(f, path) ) {
        return 1;
    }

    /// @cond
    #undef NBIPATH
    #define NBIPATH "/nbi/nbi_XXXXXXXXXX/"
    /// @endcond

    int ninj;
    if( hdf5_read_int(NBIPATH "ninj", &ninj,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Loop over injectors twice: first time we initialize offload data */
    int nbeamlet_total = 0;
    int* id = (int*) malloc( ninj*sizeof(int) );
    int* anum = (int*) malloc( ninj*sizeof(int) );
    int* znum = (int*) malloc( ninj*sizeof(int) );
    int* nbeamlet = (int*) malloc( ninj*sizeof(int) );
    real* mass = (real*) malloc( ninj*sizeof(real) );
    real* power = (real*) malloc( ninj*sizeof(real) );
    real* efrac = (real*) malloc( ninj*sizeof(real) );
    real* energy = (real*) malloc( ninj*sizeof(real) );
    real* div_h = (real*) malloc( ninj*sizeof(real) );
    real* div_v = (real*) malloc( ninj*sizeof(real) );
    real* div_halo_h = (real*) malloc( ninj*sizeof(real) );
    real* div_halo_v = (real*) malloc( ninj*sizeof(real) );
    real* div_halo_frac = (real*) malloc( ninj*sizeof(real) );
    for(int i = 0; i < ninj; i++) {
        sprintf(path, NBIPATH "inj%d/%s", i+1, "ids");
        if( hdf5_read_int(path, &id[i],
                          f, qid, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "nbeamlet");
        if( hdf5_read_int(path, &nbeamlet[i],
                          f, qid, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "power");
        if( hdf5_read_double(path, &power[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "energy");
        if( hdf5_read_double(path, &energy[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "efrac");
        if( hdf5_read_double(path, &efrac[i*3],
                             f, qid, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "divh");
        if( hdf5_read_double(path, &div_h[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divv");
        if( hdf5_read_double(path, &div_v[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divhalofrac");
        if( hdf5_read_double(path, &div_halo_frac[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divhaloh");
        if( hdf5_read_double(path, &div_halo_h[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "divhalov");
        if( hdf5_read_double(path, &div_halo_v[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}

        sprintf(path, NBIPATH "inj%d/%s", i+1, "anum");
        if( hdf5_read_int(path, &anum[i],
                          f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "znum");
        if( hdf5_read_int(path, &znum[i],
                          f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "mass");
        if( hdf5_read_double(path, &mass[i],
                             f, qid, __FILE__, __LINE__) ) {return 1;}

        /* Conver to SI */
        energy[i] *= CONST_E;
        mass[i]   *= CONST_U;

        nbeamlet_total += 6 * nbeamlet[i];
    }

    /* In second loop we now initialize the offload array which we are
     * able to allocate now */
    real* beamlet_xyz = (real*) malloc( nbeamlet_total * sizeof(real) );
    int idx = 0;
    for(int i = 0; i < data->ninj; i++) {
        int n = nbeamlet[i];
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletx");
        if( hdf5_read_double(path, &beamlet_xyz[idx + 0*n],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamlety");
        if( hdf5_read_double(path, &beamlet_xyz[idx + 1*n],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletz");
        if( hdf5_read_double(path, &beamlet_xyz[idx + 2*n],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletdx");
        if( hdf5_read_double(path, &beamlet_xyz[idx + 3*n],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletdy");
        if( hdf5_read_double(path, &beamlet_xyz[idx + 4*n],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        sprintf(path, NBIPATH "inj%d/%s", i+1, "beamletdz");
        if( hdf5_read_double(path, &beamlet_xyz[idx + 5*n],
                             f, qid, __FILE__, __LINE__) ) {return 1;}
        idx += 6*n;
    }

    /* Initialize the data */
    int err = nbi_init(data, ninj, id, anum, znum, mass, power, efrac, energy,
                       div_h, div_v, div_halo_v, div_halo_h, div_halo_frac,
                       nbeamlet, beamlet_xyz);
    return err;
}
