/**
 * @author
 * @file diag.c
 * @brief Interface for mid-simulation diagnostics
 *
 * Diagnostics are ASCOT output that is not collected at the end of the simulation
 * (like wall loads or endstate) but during the simulation. To implement a new
 * diagnostics, write functions "init", "update_gc", etc. specific to that diagnostic
 * into a new file and make sure the functions are called here when that diagnostic
 * is enabled. 
 *
 * Currently implemented diagnostics:
 * - dist_rzvpavpe: 4D (R,z,v_para,v_perp) distribution
 * - diag_orb: exact marker orbit in phase-space
 *   - orbit written at specific time intervals
 *   - Poincare plot, i.e., orbit written when marker crosses a defined plane
 *   - Last N time instances before marker is terminated where N is user-defined
 *
 * @todo Add debug distributions for time-steps etc.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ascot5.h"
#include "diag.h"
#include "diag_orb.h"
#include "dist_5D.h"
#include "dist_6D.h"
#include "particle.h"

/**
 * @brief Initializes offload array from offload data
 *
 * @param data initialized data
 * @param offload_array which is to be allocated and initialized
 */
void diag_init_offload(diag_offload_data* data, real** offload_array){
    /* Determine how long array we need and allocate it */
    int n = 0;

    if(data->orb_collect) {
    // Do nothing (for now)
    }

    if(data->dist5D_collect) {
        data->offload_dist5D_index = n;
        n += data->dist5D.n_r * data->dist5D.n_phi * data->dist5D.n_z
             * data->dist5D.n_vpara * data->dist5D.n_vperp;
    }

    if(data->dist6D_collect) {
        data->offload_dist6D_index = n;
        n += data->dist6D.n_r * data->dist6D.n_phi * data->dist6D.n_z
             * data->dist6D.n_vr * data->dist6D.n_vphi
             * data->dist6D.n_vz;
    }

    data->offload_array_length = n;
    *offload_array = malloc(n * sizeof(real));
    if(*offload_array == 0) {
        printf("Histogram memory allocation failed\n");
        exit(1);
    }

    memset(*offload_array, 0, n * sizeof(real));
}

/**
 * @brief Frees the offload data
 */
void diag_free_offload(diag_offload_data* data, real** offload_array) {
    if(data->orb_collect) {
    // Nothing to be freed really
    }
    if(data->dist5D_collect) {
    // Nothing to be freed really
    }

    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initializes diagnostics from offload data
 *
 */
void diag_init(diag_data* data, diag_offload_data* offload_data,
               real* offload_array) {
    data->diag_debug_collect = offload_data->debug_collect;
    data->diag_orb_collect = offload_data->orb_collect;
    data->diag_dist5D_collect = offload_data->dist5D_collect;
    data->diag_dist6D_collect = offload_data->dist6D_collect;

    if(data->diag_orb_collect) {
        diag_orb_init(&data->orbits, &offload_data->orbits);
    }

    if(data->diag_dist5D_collect) {
        dist_5D_init(&data->dist5D, &offload_data->dist5D,
                     &offload_array[offload_data->offload_dist5D_index]);
    }

    if(data->diag_dist6D_collect) {
        dist_6D_init(&data->dist6D, &offload_data->dist6D,
                     &offload_array[offload_data->offload_dist6D_index]);
    }
}

/**
 * @brief Collects diagnostics when marker represents a particle
 */
void diag_update_fo(diag_data* d, diag_storage* ds, particle_simd_fo* p_f,
                    particle_simd_fo* p_i) {
    if(d->diag_orb_collect) {
        diag_orb_update_fo(ds->particleId, ds->prevWriteTime, ds->nextN,
                           ds->Nlist, &d->orbits, p_f, p_i);
    }

    if(d->diag_debug_collect) {
    }

    if(d->diag_dist5D_collect) {
        dist_5D_update_fo(&d->dist5D, p_f, p_i);
    }

    if(d->diag_dist6D_collect) {
        dist_6D_update_fo(&d->dist6D, p_f, p_i);
    }
}

/**
 * @brief Collects diagnostics when marker represents a guiding center
 */
void diag_update_gc(diag_data* d, diag_storage* ds, particle_simd_gc* p_f,
                    particle_simd_gc* p_i) {
    if(d->diag_orb_collect) {
        diag_orb_update_gc(ds->particleId, ds->prevWriteTime, ds->nextN,
                           ds->Nlist, &d->orbits, p_f, p_i);
    }

    if(d->diag_debug_collect){
    }

    if(d->diag_dist5D_collect){
        dist_5D_update_gc(&d->dist5D, p_f, p_i);
    }

    if(d->diag_dist6D_collect){
        dist_6D_update_gc(&d->dist6D, p_f, p_i);
    }
}

/**
 * @brief Collects diagnostics when marker represents a magnetic field line
 */
void diag_update_ml(diag_data* d, diag_storage* ds, particle_simd_ml* p_f,
                    particle_simd_ml* p_i) {
    if(d->diag_orb_collect) {
        diag_orb_update_ml(ds->particleId, ds->prevWriteTime, ds->nextN,
                           ds->Nlist, &d->orbits, p_f, p_i);
    }

    if(d->diag_debug_collect) {
    }

    if(d->diag_dist5D_collect) {
    }
}

/**
 * @brief Sum offload data arrays as one
 * @todo This whole scheme does not work for orbits
 */
void diag_sum(diag_data* d, real* array1, real* array2) {
    if(d->diag_orb_collect) {
    // Do nothing
    }

    if(d->diag_debug_collect) {
    // Do nothing
    }

    if(d->diag_dist5D_collect){
        int start = d->offload_dist5D_index;
        int stop = start + d->dist5D.n_r * d->dist5D.n_z
                   * d->dist5D.n_vpara * d->dist5D.n_vperp;
        dist_5D_sum(start, stop, array1, array2);
    }

    if(d->diag_dist6D_collect){
        int start = d->offload_dist6D_index;
        int stop = start + d->dist6D.n_r * d->dist6D.n_phi * d->dist6D.n_z
                   * d->dist6D.n_vr * d->dist6D.n_vphi * d->dist6D.n_vz;
        dist_6D_sum(start, stop, array1, array2);
    }
}

/**
 * @brief Some diagnostics use dynamically allocated data structs which
 *        are freed here
 *
 */
void diag_clean(diag_data* d) {
    if(d->diag_orb_collect) {
        diag_orb_clean(&d->orbits);
    }

    if(d->diag_debug_collect) {
    // Do nothing
    }

    if(d->diag_dist5D_collect) {
    // Do nothing
    }
}

/**
 * @brief Allocates a new diag_storage
 */
void diag_storage_aquire(diag_data* data, diag_storage** ds) {
    *ds = malloc(sizeof(diag_storage));

    /* N-last specific input */
    (*ds)->Nlist = NULL;
    if(data->diag_orb_collect && data->orbits.mode == DIAG_ORB_WRITELAST) {
	(*ds)->Nlist = malloc(data->orbits.writeNlast*NSIMD*sizeof(diag_orb_dat*));
    }

    /* Initialize data storage */
    for(int i=0; i < NSIMD; i++) {
    (*ds)->particleId[i] = -2;
    }
}

/**
 * @brief De-allocates a diag_storage
 */
void diag_storage_discard(diag_storage* ds) {
    if(ds->Nlist != NULL) {
    free(ds->Nlist);
    }
    free(ds);
}
