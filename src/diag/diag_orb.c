/**
 * @file diag_orb.c
 * @brief Functions to write particle and guiding center information.
 */
#include "diag_orb.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ascot5.h"
#include "consts.h"
#include "options.h"
#include "diag.h"

/**
 * @brief Initializes orbit diagnostics data.
 *
 * The offload array should have a length of Nfld * Nmrk * Npnt and elements
 * initialized to zero. The orbit data will be stored in this field as
 * arr[i_field*N_markers*N_points + i_marker*N_points + i_point].
 *
 * Note that not all markers fill all space assigned to them before their
 * simulation is terminated.
 *
 * @param data orbit diagnostics data struct
 */
void diag_orb_init(diag_orb_data* data, sim_parameters* params,
                   size_t nmarkers) {

    int step = nmarkers*params->number_of_points_per_marker;
    data->id = (real*) calloc( step, sizeof(real) );
    if(params->poincare) {
        data->pncrid = (real*) calloc( step, sizeof(real) );
        data->pncrdi = (real*) calloc( step, sizeof(real) );
    }
    switch(params->simulation_mode) {

        case simulate_mode_fo:
            data->mileage = (real*) calloc( step, sizeof(real) );
            data->r       = (real*) calloc( step, sizeof(real) );
            data->phi     = (real*) calloc( step, sizeof(real) );
            data->z       = (real*) calloc( step, sizeof(real) );
            data->p_r     = (real*) calloc( step, sizeof(real) );
            data->p_phi   = (real*) calloc( step, sizeof(real) );
            data->p_z     = (real*) calloc( step, sizeof(real) );
            data->weight  = (real*) calloc( step, sizeof(real) );
            data->charge  = (real*) calloc( step, sizeof(real) );
            data->rho     = (real*) calloc( step, sizeof(real) );
            data->theta   = (real*) calloc( step, sizeof(real) );
            data->B_r     = (real*) calloc( step, sizeof(real) );
            data->B_phi   = (real*) calloc( step, sizeof(real) );
            data->B_z     = (real*) calloc( step, sizeof(real) );
            data->simmode = (real*) calloc( step, sizeof(real) );
            break;

        case simulate_mode_gc:
            data->mileage = (real*) calloc( step, sizeof(real) );
            data->r       = (real*) calloc( step, sizeof(real) );
            data->phi     = (real*) calloc( step, sizeof(real) );
            data->z       = (real*) calloc( step, sizeof(real) );
            data->ppar    = (real*) calloc( step, sizeof(real) );
            data->mu      = (real*) calloc( step, sizeof(real) );
            data->zeta    = (real*) calloc( step, sizeof(real) );
            data->weight  = (real*) calloc( step, sizeof(real) );
            data->charge  = (real*) calloc( step, sizeof(real) );
            data->rho     = (real*) calloc( step, sizeof(real) );
            data->theta   = (real*) calloc( step, sizeof(real) );
            data->B_r     = (real*) calloc( step, sizeof(real) );
            data->B_phi   = (real*) calloc( step, sizeof(real) );
            data->B_z     = (real*) calloc( step, sizeof(real) );
            data->simmode = (real*) calloc( step, sizeof(real) );
            break;

        case simulate_mode_ml:
            data->mileage = (real*) calloc( step, sizeof(real) );
            data->r       = (real*) calloc( step, sizeof(real) );
            data->phi     = (real*) calloc( step, sizeof(real) );
            data->z       = (real*) calloc( step, sizeof(real) );
            data->rho     = (real*) calloc( step, sizeof(real) );
            data->theta   = (real*) calloc( step, sizeof(real) );
            data->B_r     = (real*) calloc( step, sizeof(real) );
            data->B_phi   = (real*) calloc( step, sizeof(real) );
            data->B_z     = (real*) calloc( step, sizeof(real) );
            data->simmode = (real*) calloc( step, sizeof(real) );
            break;

        case simulate_mode_hybrid:
            data->mileage = (real*) calloc( step, sizeof(real) );
            data->r       = (real*) calloc( step, sizeof(real) );
            data->phi     = (real*) calloc( step, sizeof(real) );
            data->z       = (real*) calloc( step, sizeof(real) );
            data->p_r     = (real*) calloc( step, sizeof(real) );
            data->p_phi   = (real*) calloc( step, sizeof(real) );
            data->p_z     = (real*) calloc( step, sizeof(real) );
            data->ppar    = (real*) calloc( step, sizeof(real) );
            data->mu      = (real*) calloc( step, sizeof(real) );
            data->zeta    = (real*) calloc( step, sizeof(real) );
            data->weight  = (real*) calloc( step, sizeof(real) );
            data->charge  = (real*) calloc( step, sizeof(real) );
            data->rho     = (real*) calloc( step, sizeof(real) );
            data->theta   = (real*) calloc( step, sizeof(real) );
            data->B_r     = (real*) calloc( step, sizeof(real) );
            data->B_phi   = (real*) calloc( step, sizeof(real) );
            data->B_z     = (real*) calloc( step, sizeof(real) );
            data->simmode = (real*) calloc( step, sizeof(real) );
            break;
    }

    data->mrk_pnt = (integer*) malloc( nmarkers*sizeof(integer) );
    data->mrk_recorded = (real*) malloc( nmarkers*sizeof(real) );

    memset(data->mrk_pnt, 0, nmarkers*sizeof(integer));
    memset(data->mrk_recorded, 0, nmarkers*sizeof(real));
}

/**
 * @brief Free allocated resources
 *
 * @param data orbit diagnostics data struct
 */
void diag_orb_free(diag_orb_data* data, sim_parameters* params){
    free(data->mrk_pnt);
    free(data->mrk_recorded);

    free(data->id);
    if(params->poincare) {
        free(data->pncrid);
        free(data->pncrdi);
    }
    switch(params->record_mode) {

        case simulate_mode_fo:
            free(data->mileage);
            free(data->r);
            free(data->phi);
            free(data->z);
            free(data->p_r);
            free(data->p_phi);
            free(data->p_z);
            free(data->weight);
            free(data->charge );
            free(data->rho);
            free(data->theta);
            free(data->B_r);
            free(data->B_phi);
            free(data->B_z);
            free(data->simmode);
            break;

        case simulate_mode_gc:
            free(data->mileage);
            free(data->r);
            free(data->phi);
            free(data->z);
            free(data->ppar);
            free(data->mu);
            free(data->zeta);
            free(data->weight);
            free(data->charge );
            free(data->rho);
            free(data->theta);
            free(data->B_r);
            free(data->B_phi);
            free(data->B_z);
            free(data->simmode);
            break;

        case simulate_mode_ml:
            free(data->mileage);
            free(data->r);
            free(data->phi);
            free(data->z);
            free(data->rho);
            free(data->theta);
            free(data->B_r);
            free(data->B_phi);
            free(data->B_z);
            free(data->simmode);
            break;

        case simulate_mode_hybrid:
            free(data->mileage);
            free(data->r);
            free(data->phi);
            free(data->z);
            free(data->p_r);
            free(data->p_phi);
            free(data->p_z);
            free(data->ppar);
            free(data->mu);
            free(data->zeta);
            free(data->weight);
            free(data->charge );
            free(data->rho);
            free(data->theta);
            free(data->B_r);
            free(data->B_phi);
            free(data->B_z);
            free(data->simmode);
            break;
    }
}

/**
 * @brief Collects orbit diagnostics when marker represents a particle
 *
 * @param data orbit diagnostics data struct
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void diag_orb_update_fo(diag_orb_data* data, sim_parameters* params,
                        particle_simd_fo* p_f, particle_simd_fo* p_i) {

    if(!params->poincare) {

        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {

            /* Mask dummy markers */
            if(p_f->id[i] > 0) {

                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * params->number_of_points_per_marker
                               + ipoint;

                /* If this is the first time-step, record marker position. */
                if( data->id[imrk * params->number_of_points_per_marker] == 0 ) {
                    data->id[idx]     = (real)p_i->id[i];
                    data->mileage[idx]= p_i->mileage[i];
                    data->r[idx]      = p_i->r[i];
                    data->phi[idx]    = p_i->phi[i];
                    data->z[idx]      = p_i->z[i];
                    data->p_r[idx]    = p_i->p_r[i];
                    data->p_phi[idx]  = p_i->p_phi[i];
                    data->p_z[idx]    = p_i->p_z[i];
                    data->weight[idx] = p_i->weight[i];
                    data->charge[idx] = p_i->charge[i];
                    data->rho[idx]    = p_i->rho[i];
                    data->theta[idx]  = p_i->theta[i];
                    data->B_r[idx]    = p_i->B_r[i];
                    data->B_phi[idx]  = p_i->B_phi[i];
                    data->B_z[idx]    = p_i->B_z[i];
                    data->simmode[idx]= DIAG_ORB_FO;

                    ipoint++;
                    if(ipoint == params->number_of_points_per_marker) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_i->mileage[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + params->interval
                    - p_f->mileage[i];
                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * params->number_of_points_per_marker + ipoint;

                    data->id[idx]     = (real)p_f->id[i];
                    data->mileage[idx]= p_f->mileage[i];
                    data->r[idx]      = p_f->r[i];
                    data->phi[idx]    = p_f->phi[i];
                    data->z[idx]      = p_f->z[i];
                    data->p_r[idx]    = p_f->p_r[i];
                    data->p_phi[idx]  = p_f->p_phi[i];
                    data->p_z[idx]    = p_f->p_z[i];
                    data->weight[idx] = p_f->weight[i];
                    data->charge[idx] = p_f->charge[i];
                    data->rho[idx]    = p_f->rho[i];
                    data->theta[idx]  = p_f->theta[i];
                    data->B_r[idx]    = p_f->B_r[i];
                    data->B_phi[idx]  = p_f->B_phi[i];
                    data->B_z[idx]    = p_f->B_z[i];
                    data->simmode[idx]= DIAG_ORB_FO;

                    ipoint++;
                    if(ipoint == params->number_of_points_per_marker) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk] = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }
            }
        }
    }
    else {

        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {
            /* Mask dummy markers and those whose time-step was rejected. */
            if( p_f->id[i] > 0 && (p_f->mileage[i] != p_i->mileage[i]) ) {

                real k;
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * params->number_of_points_per_marker
                               + ipoint;

                /* Check and store toroidal crossings. */
                for(int j=0; j < params->ntoroidalplots; j++) {
                    k = diag_orb_check_plane_crossing(
                        p_f->phi[i], p_i->phi[i], params->toroidal_angles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->p_r[idx]    = k*p_f->p_r[i]    + d*p_i->p_r[i];
                        data->p_phi[idx]  = k*p_f->p_phi[i]  + d*p_i->p_phi[i];
                        data->p_z[idx]    = k*p_f->p_z[i]    + d*p_i->p_z[i];
                        data->weight[idx] = k*p_f->weight[i] + d*p_i->weight[i];
                        data->charge[idx] = p_i->charge[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] = j;
                        data->pncrdi[idx] = 1 - 2 * (p_f->phi[i] < p_i->phi[i]);
                        data->simmode[idx]= DIAG_ORB_FO;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store poloidal crossings. */
                for(int j=0; j < params->npoloidalplots; j++) {
                    k = diag_orb_check_plane_crossing(
                        p_f->theta[i], p_i->theta[i], params->poloidal_angles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->p_r[idx]    = k*p_f->p_r[i]    + d*p_i->p_r[i];
                        data->p_phi[idx]  = k*p_f->p_phi[i]  + d*p_i->p_phi[i];
                        data->p_z[idx]    = k*p_f->p_z[i]    + d*p_i->p_z[i];
                        data->weight[idx] = k*p_f->weight[i] + d*p_i->weight[i];
                        data->charge[idx] = p_i->charge[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] = j + params->ntoroidalplots;
                        data->pncrdi[idx] = 1 - 2 * (p_f->theta[i] < p_i->theta[i]);
                        data->simmode[idx]= DIAG_ORB_FO;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store radial crossings. */
                for(int j=0; j < params->nradialplots; j++) {
                    k = diag_orb_check_radial_crossing(
                        p_f->rho[i],p_i->rho[i], params->radial_distances[j]);
                    if(k) {
                        real d = k;
                        k = 1-d;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->p_r[idx]    = k*p_f->p_r[i]    + d*p_i->p_r[i];
                        data->p_phi[idx]  = k*p_f->p_phi[i]  + d*p_i->p_phi[i];
                        data->p_z[idx]    = k*p_f->p_z[i]    + d*p_i->p_z[i];
                        data->weight[idx] = k*p_f->weight[i] + d*p_i->weight[i];
                        data->charge[idx] = p_i->charge[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] = j + params->ntoroidalplots + params->npoloidalplots;
                        data->pncrdi[idx] = 1 - 2 * (p_f->rho[i] < p_i->rho[i]);
                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }
            }
        }
    }
}

/**
 * @brief Collects orbit diagnostics when marker represents a guiding center
 *
 * @param data orbit diagnostics data struct
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void diag_orb_update_gc(diag_orb_data* data, sim_parameters* params,
                        particle_simd_gc* p_f, particle_simd_gc* p_i) {

    if(!params->poincare) {
        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {

            /* Mask dummy markers */
            if(p_f->id[i] > 0) {
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * params->number_of_points_per_marker + ipoint;

                /* If this is the first time-step, record marker position. */
                if( data->id[imrk * params->number_of_points_per_marker] == 0 ) {
                    data->id[idx]     = (real)(p_i->id[i]);
                    data->mileage[idx]= p_i->mileage[i];
                    data->r[idx]      = p_i->r[i];
                    data->phi[idx]    = p_i->phi[i];
                    data->z[idx]      = p_i->z[i];
                    data->ppar[idx]   = p_i->ppar[i];
                    data->mu[idx]     = p_i->mu[i];
                    data->zeta[idx]   = p_i->zeta[i];
                    data->weight[idx] = p_i->weight[i];
                    data->charge[idx] = p_i->charge[i];
                    data->rho[idx]    = p_i->rho[i];
                    data->theta[idx]  = p_i->theta[i];
                    data->B_r[idx]    = p_i->B_r[i];
                    data->B_phi[idx]  = p_i->B_phi[i];
                    data->B_z[idx]    = p_i->B_z[i];
                    data->simmode[idx]= DIAG_ORB_GC;

                    ipoint++;
                    if(ipoint == params->number_of_points_per_marker) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_i->mileage[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + params->interval
                    - p_f->mileage[i];

                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * params->number_of_points_per_marker + ipoint;

                    data->id[idx]     = (real)p_f->id[i];
                    data->mileage[idx]= p_f->mileage[i];
                    data->r[idx]      = p_f->r[i];
                    data->phi[idx]    = p_f->phi[i];
                    data->z[idx]      = p_f->z[i];
                    data->ppar[idx]   = p_f->ppar[i];
                    data->mu[idx]     = p_f->mu[i];
                    data->zeta[idx]   = p_f->zeta[i];
                    data->weight[idx] = p_f->weight[i];
                    data->charge[idx] = p_f->charge[i];
                    data->rho[idx]    = p_f->rho[i];
                    data->theta[idx]  = p_f->theta[i];
                    data->B_r[idx]    = p_f->B_r[i];
                    data->B_phi[idx]  = p_f->B_phi[i];
                    data->B_z[idx]    = p_f->B_z[i];
                    data->simmode[idx]= DIAG_ORB_GC;

                    ipoint++;
                    if(ipoint == params->number_of_points_per_marker) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }
            }
        }
    }
    else {
        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {
            /* Mask dummy markers and those whose time-step was rejected. */
            if( p_f->id[i] > 0 && (p_f->mileage[i] != p_i->mileage[i]) ) {

                real k;
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * params->number_of_points_per_marker + ipoint;

                /* Check and store toroidal crossings. */
                for(int j=0; j < params->ntoroidalplots; j++) {
                    k = diag_orb_check_plane_crossing(
                        p_f->phi[i], p_i->phi[i], params->toroidal_angles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->ppar[idx]   = k*p_f->ppar[i]   + d*p_i->ppar[i];
                        data->mu[idx]     = k*p_f->mu[i]     + d*p_i->mu[i];
                        data->zeta[idx]   = k*p_f->zeta[i]   + d*p_i->zeta[i];
                        data->weight[idx] = k*p_f->weight[i] + d*p_i->weight[i];
                        data->charge[idx] = p_i->charge[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] = j;
                        data->pncrdi[idx] = 1 - 2 * (p_f->phi[i] < p_i->phi[i]);
                        data->simmode[idx]= DIAG_ORB_GC;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store poloidal crossings. */
                for(int j=0; j < params->npoloidalplots; j++) {
                    k = diag_orb_check_plane_crossing(
                        p_f->theta[i], p_i->theta[i], params->poloidal_angles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->ppar[idx]   = k*p_f->ppar[i]   + d*p_i->ppar[i];
                        data->mu[idx]     = k*p_f->mu[i]     + d*p_i->mu[i];
                        data->zeta[idx]   = k*p_f->zeta[i]   + d*p_i->zeta[i];
                        data->weight[idx] = k*p_f->weight[i] + d*p_i->weight[i];
                        data->charge[idx] = p_i->charge[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] = j + params->ntoroidalplots;
                        data->pncrdi[idx] = 1 - 2 * (p_f->theta[i] < p_i->theta[i]);
                        data->simmode[idx]= DIAG_ORB_GC;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }


                /* Check and store radial crossings. */
                for(int j=0; j < params->nradialplots; j++) {
                    k = diag_orb_check_radial_crossing(
                        p_f->rho[i], p_i->rho[i], params->radial_distances[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->ppar[idx]   = k*p_f->ppar[i]   + d*p_i->ppar[i];
                        data->mu[idx]     = k*p_f->mu[i]     + d*p_i->mu[i];
                        data->zeta[idx]   = k*p_f->zeta[i]   + d*p_i->zeta[i];
                        data->weight[idx] = k*p_f->weight[i] + d*p_i->weight[i];
                        data->charge[idx] = p_i->charge[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] =
                                j + params->ntoroidalplots + params->npoloidalplots;
                        data->pncrdi[idx] = 1 - 2 * (p_f->rho[i] < p_i->rho[i]);

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }
            }
        }
    }
}

/**
 * @brief Collects orbit diagnostics when marker represents a field line
 *
 * @param data orbit diagnostics data struct
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void diag_orb_update_ml(diag_orb_data* data, sim_parameters* params,
                        particle_simd_ml* p_f, particle_simd_ml* p_i) {

    if(!params->poincare) {

        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {

            /* Mask dummy markers */
            if(p_f->id[i] > 0) {
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * params->number_of_points_per_marker + ipoint;

                /* If this is the first time-step, record marker position. */
                if( data->id[imrk * params->number_of_points_per_marker] == 0 ) {
                    data->id[idx]      = (real)p_i->id[i];
                    data->mileage[idx] = p_i->mileage[i];
                    data->r[idx]       = p_i->r[i];
                    data->phi[idx]     = p_i->phi[i];
                    data->z[idx]       = p_i->z[i];
                    data->rho[idx]     = p_i->rho[i];
                    data->theta[idx]   = p_i->theta[i];
                    data->B_r[idx]     = p_i->B_r[i];
                    data->B_phi[idx]   = p_i->B_phi[i];
                    data->B_z[idx]     = p_i->B_z[i];
                    data->simmode[idx] = DIAG_ORB_ML;

                    ipoint++;
                    if(ipoint == params->number_of_points_per_marker) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_i->mileage[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + params->interval
                    - p_f->mileage[i];
                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * params->number_of_points_per_marker + ipoint;
                    data->id[idx]      = (real)p_f->id[i];
                    data->mileage[idx] = p_f->mileage[i];
                    data->r[idx]       = p_f->r[i];
                    data->phi[idx]     = p_f->phi[i];
                    data->z[idx]       = p_f->z[i];
                    data->rho[idx]     = p_f->rho[i];
                    data->theta[idx]   = p_f->theta[i];
                    data->B_r[idx]     = p_f->B_r[i];
                    data->B_phi[idx]   = p_f->B_phi[i];
                    data->B_z[idx]     = p_f->B_z[i];
                    data->simmode[idx] = DIAG_ORB_ML;

                    ipoint++;
                    if(ipoint == params->number_of_points_per_marker) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }
            }
        }
    }
    else {
        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {
            /* Mask dummy markers and thosw whose time-step was rejected. */
            if( p_f->id[i] > 0 && (p_f->mileage[i] != p_i->mileage[i]) ) {

                real k;
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * params->number_of_points_per_marker + ipoint;

                /* Check and store toroidal crossings. */
                for(int j=0; j < params->ntoroidalplots; j++) {
                    k = diag_orb_check_plane_crossing(
                        p_f->phi[i], p_i->phi[i], params->toroidal_angles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i]+ d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]      + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]    + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]      + d*p_i->z[i];
                        data->rho[idx]    = k*p_f->rho[i]    + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]  + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]    + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]  + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]    + d*p_i->B_z[i];
                        data->pncrid[idx] = j;
                        data->pncrdi[idx] = 1 - 2 * (p_f->phi[i] < p_i->phi[i]);
                        data->simmode[idx]= DIAG_ORB_ML;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store poloidal crossings. */
                for(int j=0; j < params->npoloidalplots; j++) {
                    k = diag_orb_check_plane_crossing(
                        p_f->theta[i], p_i->theta[i], params->poloidal_angles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i] + d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]       + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]     + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]       + d*p_i->z[i];
                        data->rho[idx]    = k*p_f->rho[i]     + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]   + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]     + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]   + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]     + d*p_i->B_z[i];
                        data->pncrid[idx] = j + params->ntoroidalplots;
                        data->pncrdi[idx] = 1 - 2 * (p_f->theta[i] < p_i->theta[i]);
                        data->simmode[idx]= DIAG_ORB_ML;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store radial crossings. */
                for(int j=0; j < params->nradialplots; j++) {
                    k = diag_orb_check_radial_crossing(
                        p_f->rho[i], p_i->rho[i], params->radial_distances[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * params->number_of_points_per_marker + ipoint;
                        data->id[idx]     = (real)p_f->id[i];
                        data->mileage[idx]= k*p_f->mileage[i] + d*p_i->mileage[i];
                        data->r[idx]      = k*p_f->r[i]       + d*p_i->r[i];
                        data->phi[idx]    = k*p_f->phi[i]     + d*p_i->phi[i];
                        data->z[idx]      = k*p_f->z[i]       + d*p_i->z[i];
                        data->rho[idx]    = k*p_f->rho[i]     + d*p_i->rho[i];
                        data->theta[idx]  = k*p_f->theta[i]   + d*p_i->theta[i];
                        data->B_r[idx]    = k*p_f->B_r[i]     + d*p_i->B_r[i];
                        data->B_phi[idx]  = k*p_f->B_phi[i]   + d*p_i->B_phi[i];
                        data->B_z[idx]    = k*p_f->B_z[i]     + d*p_i->B_z[i];
                        data->pncrid[idx] =
                                j + params->ntoroidalplots + params->npoloidalplots;

                        ipoint++;
                        if(ipoint == params->number_of_points_per_marker) {
                            ipoint = 0;
                        }

                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }
            }
        }
    }
}

/**
 * @brief Check if marker has crossed a plane.
 *
 * This helper function checks whether the angle, either toroidal or poloidal,
 * that defines a Poincare plane is between marker's initial and final angles
 * (of single timestep).
 *
 * @param fang marker final angle in radians.
 * @param iang marker initial angle in radians.
 * @param ang0 Poincare plane angle.
 *
 * @return zero if no-crossing, number k, ang0 = k + (fang - iang), otherwise.
 */
real diag_orb_check_plane_crossing(real fang, real iang, real ang0) {

    real k = 0;
    /* Check whether nag0 is between iang and fang. Note that this     *
     * implementation works only because iang and fang are cumulative. */
    if( floor( (fang - ang0)/CONST_2PI ) != floor( (iang - ang0)/CONST_2PI ) ) {

        /* Move iang to interval [0, 2pi] */
        real a = fmod(iang, CONST_2PI);
        if(a < 0){
            a = CONST_2PI + a;
        }

        a = fabs(ang0 - a);
        if(a > CONST_PI) {
            a = CONST_2PI - a;
        }
        k = fabs(a / (fang - iang));
    }
    return k;
}

/**
 * @brief Check if marker has crossed given rho
 *
 * This helper function checks whether given rho that defines a Poincare plane
 * is between marker's initial and final rhos (of single timestep).
 *
 * @param frho marker final rho in metres.
 * @param irho marker initial rho in metres.
 * @param rho0 Poincare plane rho.
 *
 * @return zero if no-crossing, number k, rho0 = k * (frho - irho), otherwise.
 */
real diag_orb_check_radial_crossing(real frho, real irho, real rho0){

    real k = 0;
    if((frho <= rho0 && irho > rho0) || (irho <= rho0 && frho > rho0)){
        k = fabs((irho - rho0) / (frho - irho));
    }
    return k;
}
