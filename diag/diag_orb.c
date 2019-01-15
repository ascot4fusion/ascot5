/**
 * @file diag_orb.c
 * @brief Functions to write particle and guiding center information.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../simulate.h"
#include "../gctransform.h"
#include "../particle.h"
#include "../B_field.h"
#include "diag_orb.h"

/**
 * @brief Initializes orbit diagnostics offload data.
 *
 * The offload array should have a length of Nfld * Nmrk * Npnt and elements
 * initialized to zero. The orbit data will be stored in this field as
 * arr[i_field*N_markers*N_points + i_marker*N_points + i_point].
 *
 * Note that not all markers fill all space assigned to them before their
 * simulation is terminated.
 *
 * @param data orbit diagnostics data struct
 * @param offload_data orbit diagnostics offload data struct
 * @param offload_array offload data array
 */
void diag_orb_init(diag_orb_data* data, diag_orb_offload_data* offload_data,
                   real* offload_array) {

    data->mode = offload_data->mode;
    data->Nmrk = offload_data->Nmrk;
    data->Npnt = offload_data->Npnt;

    int step = data->Nmrk*data->Npnt;

    if(data->mode == DIAG_ORB_INTERVAL) {
        data->writeInterval = offload_data->writeInterval;
    }
    else if(data->mode == DIAG_ORB_POINCARE) {
        data->ntoroidalplots = offload_data->ntoroidalplots;
        for(int i=0; i<data->ntoroidalplots; i++) {
            data->toroidalangles[i] = offload_data->toroidalangles[i];
        }
        data->npoloidalplots = offload_data->npoloidalplots;
        for(int i=0; i<data->ntoroidalplots; i++) {
            data->poloidalangles[i] = offload_data->poloidalangles[i];
        }

        data->pncrid = &(offload_array[step*offload_data->Nfld]);
    }

    switch(offload_data->record_mode) {

        case simulate_mode_fo:
            data->id     = &(offload_array[step*0]);
            data->time   = &(offload_array[step*1]);
            data->r      = &(offload_array[step*2]);
            data->phi    = &(offload_array[step*3]);
            data->z      = &(offload_array[step*4]);
            data->rdot   = &(offload_array[step*5]);
            data->phidot = &(offload_array[step*6]);
            data->zdot   = &(offload_array[step*7]);
            data->weight = &(offload_array[step*8]);
            data->charge = &(offload_array[step*9]);
            data->rho    = &(offload_array[step*10]);
            data->pol    = &(offload_array[step*11]);
            data->B_r    = &(offload_array[step*12]);
            data->B_phi  = &(offload_array[step*13]);
            data->B_z    = &(offload_array[step*14]);
            break;

        case simulate_mode_gc:
            data->id     = &(offload_array[step*0]);
            data->time   = &(offload_array[step*1]);
            data->r      = &(offload_array[step*2]);
            data->phi    = &(offload_array[step*3]);
            data->z      = &(offload_array[step*4]);
            data->vpar   = &(offload_array[step*5]);
            data->mu     = &(offload_array[step*6]);
            data->theta  = &(offload_array[step*7]);
            data->weight = &(offload_array[step*8]);
            data->charge = &(offload_array[step*9]);
            data->rho    = &(offload_array[step*10]);
            data->pol    = &(offload_array[step*11]);
            data->B_r    = &(offload_array[step*12]);
            data->B_phi  = &(offload_array[step*13]);
            data->B_z    = &(offload_array[step*14]);
            break;

        case simulate_mode_ml:
            data->id     = &(offload_array[step*0]);
            data->time   = &(offload_array[step*1]);
            data->r      = &(offload_array[step*2]);
            data->phi    = &(offload_array[step*3]);
            data->z      = &(offload_array[step*4]);
            data->rho    = &(offload_array[step*5]);
            data->pol    = &(offload_array[step*6]);
            data->B_r    = &(offload_array[step*7]);
            data->B_phi  = &(offload_array[step*8]);
            data->B_z    = &(offload_array[step*9]);
            break;

        case simulate_mode_hybrid:
            data->id     = &(offload_array[step*0]);
            data->time   = &(offload_array[step*1]);
            data->r      = &(offload_array[step*2]);
            data->phi    = &(offload_array[step*3]);
            data->z      = &(offload_array[step*4]);
            data->rdot   = &(offload_array[step*5]);
            data->phidot = &(offload_array[step*6]);
            data->zdot   = &(offload_array[step*7]);
            data->vpar   = &(offload_array[step*8]);
            data->mu     = &(offload_array[step*9]);
            data->theta  = &(offload_array[step*10]);
            data->weight = &(offload_array[step*11]);
            data->charge = &(offload_array[step*12]);
            data->rho    = &(offload_array[step*13]);
            data->pol    = &(offload_array[step*14]);
            data->B_r    = &(offload_array[step*15]);
            data->B_phi  = &(offload_array[step*16]);
            data->B_z    = &(offload_array[step*17]);
            break;
    }

    data->mrk_pnt      = malloc( data->Nmrk*sizeof(integer) );
    data->mrk_recorded = malloc( data->Nmrk*sizeof(real) );

    memset(data->mrk_pnt, 0, data->Nmrk*sizeof(integer));
    memset(data->mrk_recorded, 0, data->Nmrk*sizeof(real));
}

/**
 * @brief Free orbit diagnostics data
 *
 * @param data orbit diagnostics data struct
 */
void diag_orb_free(diag_orb_data* data){
    free(data->mrk_pnt);
    free(data->mrk_recorded);
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
void diag_orb_update_fo(diag_orb_data* data, particle_simd_fo* p_f,
                        particle_simd_fo* p_i) {

    if(data->mode == DIAG_ORB_INTERVAL) {

        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {

            /* Mask dummy markers */
            if(p_f->id[i] > 0) {
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * data->Npnt + ipoint;

                /* If this is the first time-step, record marker position. */
                if( data->id[imrk * data->Npnt] == 0 ) {
                    data->id[idx]     = (real)p_i->id[i];
                    data->time[idx]   = p_i->time[i];
                    data->r[idx]      = p_i->r[i];
                    data->phi[idx]    = p_i->phi[i];
                    data->z[idx]      = p_i->z[i];
                    data->rdot[idx]   = p_i->rdot[i];
                    data->phidot[idx] = p_i->phidot[i];
                    data->zdot[idx]   = p_i->zdot[i];
                    data->weight[idx] = p_i->weight[i];
                    data->charge[idx] = p_i->charge[i];
                    data->rho[idx]    = p_i->rho[i];
                    data->pol[idx]    = p_i->pol[i];
                    data->B_r[idx]    = p_i->B_r[i];
                    data->B_phi[idx]  = p_i->B_phi[i];
                    data->B_z[idx]    = p_i->B_z[i];

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->time[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + data->writeInterval
                    - p_f->time[i];
                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * data->Npnt + ipoint;

                    data->id[idx]     = (real)p_f->id[i];
                    data->time[idx]   = p_f->time[i];
                    data->r[idx]      = p_f->r[i];
                    data->phi[idx]    = p_f->phi[i];
                    data->z[idx]      = p_f->z[i];
                    data->rdot[idx]   = p_f->rdot[i];
                    data->phidot[idx] = p_f->phidot[i];
                    data->zdot[idx]   = p_f->zdot[i];
                    data->weight[idx] = p_f->weight[i];
                    data->charge[idx] = p_f->charge[i];
                    data->rho[idx]    = p_f->rho[i];
                    data->pol[idx]    = p_f->pol[i];
                    data->B_r[idx]    = p_f->B_r[i];
                    data->B_phi[idx]  = p_f->B_phi[i];
                    data->B_z[idx]    = p_f->B_z[i];

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->time[i];
                }
            }
        }
    }
    else if(data->mode == DIAG_ORB_POINCARE) {

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
void diag_orb_update_gc(diag_orb_data* data, particle_simd_gc* p_f,
                        particle_simd_gc* p_i) {

    if(data->mode == DIAG_ORB_INTERVAL) {
        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {

            /* Mask dummy markers */
            if(p_f->id[i] > 0) {
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * data->Npnt + ipoint;

                /* If this is the first time-step, record marker position. */
                if( data->id[imrk * data->Npnt] == 0 ) {
                    data->id[idx]     = (real)(p_i->id[i]);
                    data->time[idx]   = p_i->time[i];
                    data->r[idx]      = p_i->r[i];
                    data->phi[idx]    = p_i->phi[i];
                    data->z[idx]      = p_i->z[i];
                    data->vpar[idx]   = p_i->vpar[i];
                    data->mu[idx]     = p_i->mu[i];
                    data->theta[idx]  = p_i->theta[i];
                    data->weight[idx] = p_i->weight[i];
                    data->charge[idx] = p_i->charge[i];
                    data->rho[idx]    = p_i->rho[i];
                    data->pol[idx]    = p_i->pol[i];
                    data->B_r[idx]    = p_i->B_r[i];
                    data->B_phi[idx]  = p_i->B_phi[i];
                    data->B_z[idx]    = p_i->B_z[i];

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->time[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + data->writeInterval
                    - p_f->time[i];

                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                                    idx = imrk * data->Npnt + ipoint;

                    data->id[idx]     = (real)p_f->id[i];
                    data->time[idx]   = p_f->time[i];
                    data->r[idx]      = p_f->r[i];
                    data->phi[idx]    = p_f->phi[i];
                    data->z[idx]      = p_f->z[i];
                    data->vpar[idx]   = p_f->vpar[i];
                    data->mu[idx]     = p_f->mu[i];
                    data->theta[idx]  = p_f->theta[i];
                    data->weight[idx] = p_f->weight[i];
                    data->charge[idx] = p_f->charge[i];
                    data->rho[idx]    = p_f->rho[i];
                    data->pol[idx]    = p_f->pol[i];
                    data->B_r[idx]    = p_f->B_r[i];
                    data->B_phi[idx]  = p_f->B_phi[i];
                    data->B_z[idx]    = p_f->B_z[i];
                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->time[i];
                }
            }
        }
    }
    else if(data->mode == DIAG_ORB_POINCARE) {

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
void diag_orb_update_ml(diag_orb_data* data, particle_simd_ml* p_f,
                        particle_simd_ml* p_i) {

    if(data->mode == DIAG_ORB_INTERVAL) {

        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {

            /* Mask dummy markers */
            if(p_f->id[i] > 0) {
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * data->Npnt + ipoint;

                /* If this is the first time-step, record marker position. */
                if( data->id[imrk * data->Npnt] == 0 ) {
                    data->id[idx]    = (real)p_i->id[i];
                    data->time[idx]  = p_i->time[i];
                    data->r[idx]     = p_i->r[i];
                    data->phi[idx]   = p_i->phi[i];
                    data->z[idx]     = p_i->z[i];
                    data->rho[idx]   = p_i->rho[i];
                    data->pol[idx]   = p_i->pol[i];
                    data->B_r[idx]   = p_i->B_r[i];
                    data->B_phi[idx] = p_i->B_phi[i];
                    data->B_z[idx]   = p_i->B_z[i];

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->time[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + data->writeInterval
                    - p_f->time[i];
                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * data->Npnt + ipoint;
                    data->id[idx]    = (real)p_f->id[i];
                    data->time[idx]  = p_f->time[i];
                    data->r[idx]     = p_f->r[i];
                    data->phi[idx]   = p_f->phi[i];
                    data->z[idx]     = p_f->z[i];
                    data->rho[idx]   = p_f->rho[i];
                    data->pol[idx]   = p_f->pol[i];
                    data->B_r[idx]   = p_f->B_r[i];
                    data->B_phi[idx] = p_f->B_phi[i];
                    data->B_z[idx]   = p_f->B_z[i];

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->time[i];
                }
            }
        }
    }
    else if(data->mode == DIAG_ORB_POINCARE) {
        
    }
}

/**
 * @brief Check whether a marker qualifies for Poincare-mode writing
 *
 * The marker qualifies if following conditions are met:
 * - It is not a dummy marker
 * - The time step was accepted
 * - Marker has crossed one of the specified poloidal or toroidal planes
 *
 * Marker is assumed to cross only maximum of one poloidal and one toroidal
 * plane. Other crosses are ignored.
 * The phi coordinate we use is "unmodulated", i.e., it is not limited to interval [0,2pi).
 * We can then find whether this poloidal plane was crossed by adding that plane's toroidal
 * coordinate on marker initial and final position, and see if the division with 2pi gives the
 * same value (no crossing) or not (marker has crossed the plane)
 *
 * Checks are done for NSIMD markers simultaneously.
 */
void diag_orb_poincareTrigger(diag_orb_data* data, int* pol, int* tor,
                              real* ftime, real* fpol, real* fphi,
                              real* itime, real* ipol, real* iphi,
                              real* kpol, real* ktor, integer* id){
    #pragma omp simd
    for(int i= 0; i < NSIMD; i++) {
        pol[i] = -1;
        tor[i] = -1;

        /* Dummy markers and those whose time-step was rejected are not
           accepted*/
        if( (id[i] != -1) && (ftime[i] != itime[i]) ) {

            /* Check if the particle has crossed one of the poloidal planes */
            for(int ip = 0; ip < data->npoloidalplots; ip++) {

                if( floor( (fphi[i] + data->poloidalangles[ip])/CONST_2PI ) !=
                    floor( (iphi[i] + data->poloidalangles[ip])/CONST_2PI )
                    ) {
                    pol[i] = ip;

                    /* Angles to interval [0, 2pi] */
                    real a = fmod(iphi[i], CONST_2PI);
                    if(a < 0){a = CONST_2PI + a;}

                    a = fabs(data->poloidalangles[ip] - a);
                    if(a > CONST_PI){a = CONST_2PI - a;}
                    kpol[i] = fabs(a / (fphi[i] - iphi[i]));
                    break;
                }
            }

            /* Check if the particle has crossed one of the toroidal planes */
            for(int ip = 0; ip < data->ntoroidalplots; ip++) {

                if( floor( (fpol[i] + data->toroidalangles[ip])/CONST_2PI ) !=
                    floor( (ipol[i] + data->toroidalangles[ip])/CONST_2PI )
                    ) {
                    tor[i] = ip;

                    /* Angles to interval [0, 2pi] */
                    real a = fmod(ipol[i], CONST_2PI);
                    if(a < 0){a = CONST_2PI + a;}

                    a = fabs(data->toroidalangles[ip] - a);
                    if(a > CONST_PI){a = CONST_2PI -a;}
                    ktor[i] = fabs(a / (fpol[i] - ipol[i]));
                    break;
                }
            }
        }
    }


}
