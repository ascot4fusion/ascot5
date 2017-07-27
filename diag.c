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
#include "ascot5.h"
#include "particle.h"
#include "distributions.h"
#include "diag_orb.h"
#include "diag.h"

/**
 * @brief Initializes offload array from offload data
 *
 * @param data initialized data
 * @param offload_array which is to be allocated and initialized
 */
void diag_init_offload(diag_offload_data* data, real** offload_array){
    /* Determine how long array we need and allocate it */
    int nlen = 1;
    int nlen_dist4D;

    if(data->orb_collect) {
	// Do nothing (for now)
    }
    if(data->dist4D_collect) {
	nlen_dist4D = data->dist4D.n_r * data->dist4D.n_z
	    * data->dist4D.n_vpara * data->dist4D.n_vperp;
    }

    
    nlen = nlen + nlen_dist4D;

    data->offload_array_length = nlen;
    *offload_array = malloc(nlen*sizeof(real));
    if(*offload_array == 0) {
        printf("Histogram memory allocation failed\n");
        exit(1);
    }

    /* Initialize array */
    nlen = 0;
    
    if(data->dist4D_collect) {
	data->offload_dist4D_index = nlen;
	for(int i = nlen; i < nlen + nlen_dist4D; i++) {
	    (*offload_array)[i] = 0.0;
	}
	nlen = nlen + nlen_dist4D;
    }
}

/**
 * @brief Frees the offload data
 */
void diag_free_offload(diag_offload_data* data, real** offload_array) {
    if(data->orb_collect) {
	// Nothing to be freed really
    }
    if(data->dist4D_collect) {
	// Nothing to be freed really
    }

    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initializes diagnostics from offload data
 *
 */
void diag_init(diag_data* data, diag_offload_data* offload_data, real* offload_array){
    data->diag_debug_collect = offload_data->debug_collect;
    data->diag_orb_collect = offload_data->orb_collect;
    data->diag_dist4D_collect = offload_data->dist4D_collect;

    if(data->diag_orb_collect) {
	diag_orb_init(&data->orbits, &offload_data->orbits);
    }
    if(data->diag_dist4D_collect) {
	dist_rzvv_init(&data->dist4D, &offload_data->dist4D, 
		       &offload_array[offload_data->offload_dist4D_index]);	
    }
}

/**
 * @brief Collects diagnostics when marker represents a particle
 */
void diag_update_fo(diag_data* d, diag_storage* ds, particle_simd_fo* p_f, particle_simd_fo* p_i){
    if(d->diag_orb_collect){
	diag_orb_update_fo(ds->particleId, ds->prevWriteTime, ds->nextN, ds->Nlist, 
			   &d->orbits, p_f, p_i);
    }
    if(d->diag_debug_collect){
	
    }
    if(d->diag_dist4D_collect){
	dist_rzvv_update_fo(&d->dist4D, p_f, p_i);
    }
}

/**
 * @brief Collects diagnostics when marker represents a guiding center
 */
void diag_update_gc(diag_data* d, diag_storage* ds, particle_simd_gc* p_f, particle_simd_gc* p_i){
    
    if(d->diag_orb_collect){
	diag_orb_update_gc(ds->particleId, ds->prevWriteTime, ds->nextN, ds->Nlist, 
			   &d->orbits, p_f, p_i);
    }
    if(d->diag_debug_collect){
	
	
    }
    if(d->diag_dist4D_collect){
	dist_rzvv_update_gc(&d->dist4D, p_f, p_i);
    }
}

/**
 * @brief Collects diagnostics when marker represents a magnetic field line
 */
void diag_update_ml(diag_data* d, diag_storage* ds, particle_simd_ml* p_f, particle_simd_ml* p_i){
    if(d->diag_orb_collect){
	diag_orb_update_ml(ds->particleId, ds->prevWriteTime, ds->nextN, ds->Nlist, 
			   &d->orbits, p_f, p_i);
    }
    if(d->diag_debug_collect){
	
    }
    if(d->diag_dist4D_collect){
	
    }
}

/**
 * @brief Sum offload data arrays as one
 * @todo This whole scheme does not work for orbits
 */
void diag_sum(diag_data* d, real* array1, real* array2){
    if(d->diag_orb_collect){
	// Do nothing
    }
    if(d->diag_debug_collect){
	// Do nothing
    }
    if(d->diag_dist4D_collect){
	int start = d->offload_dist4D_index;
	int stop = start + d->dist4D.n_r * d->dist4D.n_z
	    * d->dist4D.n_vpara * d->dist4D.n_vperp;
	dist_rzvv_sum(start, stop, array1, array2);
    }
}

/**
 * @brief Some diagnostics use dynamically allocated data structs which 
 *        are freed here
 *
 */
void diag_clean(diag_data* d){

    if(d->diag_orb_collect){
	diag_orb_clean(&d->orbits);
    }
    if(d->diag_debug_collect){
	// Do nothing
    }
    if(d->diag_dist4D_collect){
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
    if(data->orbits.mode == DIAG_ORB_WRITELAST) {
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
