#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../physlib.h"
#include "../simulate.h"
#include "../gctransform.h"
#include "../particle.h"
#include "../B_field.h"
#include "diag_orb.h"

#pragma omp declare target
#pragma omp declare simd uniform(ang0)
real diag_orb_check_plane_crossing(real fang, real iang, real ang0);
real diag_orb_check_radial_crossing(real fr, real ir, real r0);
#pragma omp end declare target

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
        for(int i=0; i<data->npoloidalplots; i++) {
            data->poloidalangles[i] = offload_data->poloidalangles[i];
        }
        data->nradialplots = offload_data->nradialplots;
        for(int i=0; i<data->nradialplots; i++) {
            data->radialdistances[i] = offload_data->radialdistances[i];
        }

        data->pncrid = &(offload_array[step*offload_data->Nfld]);
    }

    switch(offload_data->record_mode) {

        case simulate_mode_fo:
            data->id      = &(offload_array[step*0]);
            data->mileage = &(offload_array[step*1]);
            data->r       = &(offload_array[step*2]);
            data->phi     = &(offload_array[step*3]);
            data->z       = &(offload_array[step*4]);
            data->p_r     = &(offload_array[step*5]);
            data->p_phi   = &(offload_array[step*6]);
            data->p_z     = &(offload_array[step*7]);
            data->weight  = &(offload_array[step*8]);
            data->charge  = &(offload_array[step*9]);
            data->rho     = &(offload_array[step*10]);
            data->theta   = &(offload_array[step*11]);
            data->B_r     = &(offload_array[step*12]);
            data->B_phi   = &(offload_array[step*13]);
            data->B_z     = &(offload_array[step*14]);
            break;

        case simulate_mode_gc:
            data->id      = &(offload_array[step*0]);
            data->mileage = &(offload_array[step*1]);
            data->r       = &(offload_array[step*2]);
            data->phi     = &(offload_array[step*3]);
            data->z       = &(offload_array[step*4]);
            data->ppar    = &(offload_array[step*5]);
            data->mu      = &(offload_array[step*6]);
            data->zeta    = &(offload_array[step*7]);
            data->weight  = &(offload_array[step*8]);
            data->charge  = &(offload_array[step*9]);
            data->rho     = &(offload_array[step*10]);
            data->theta   = &(offload_array[step*11]);
            data->B_r     = &(offload_array[step*12]);
            data->B_phi   = &(offload_array[step*13]);
            data->B_z     = &(offload_array[step*14]);
            break;

        case simulate_mode_ml:
            data->id      = &(offload_array[step*0]);
            data->mileage = &(offload_array[step*1]);
            data->r       = &(offload_array[step*2]);
            data->phi     = &(offload_array[step*3]);
            data->z       = &(offload_array[step*4]);
            data->rho     = &(offload_array[step*5]);
            data->theta   = &(offload_array[step*6]);
            data->B_r     = &(offload_array[step*7]);
            data->B_phi   = &(offload_array[step*8]);
            data->B_z     = &(offload_array[step*9]);
            break;

        case simulate_mode_hybrid:
            data->id      = &(offload_array[step*0]);
            data->mileage = &(offload_array[step*1]);
            data->r       = &(offload_array[step*2]);
            data->phi     = &(offload_array[step*3]);
            data->z       = &(offload_array[step*4]);
            data->p_r     = &(offload_array[step*5]);
            data->p_phi   = &(offload_array[step*6]);
            data->p_z     = &(offload_array[step*7]);
            data->ppar    = &(offload_array[step*8]);
            data->mu      = &(offload_array[step*9]);
            data->zeta    = &(offload_array[step*10]);
            data->weight  = &(offload_array[step*11]);
            data->charge  = &(offload_array[step*12]);
            data->rho     = &(offload_array[step*13]);
            data->theta   = &(offload_array[step*14]);
            data->B_r     = &(offload_array[step*15]);
            data->B_phi   = &(offload_array[step*16]);
            data->B_z     = &(offload_array[step*17]);
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

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + data->writeInterval
                    - p_f->mileage[i];
                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * data->Npnt + ipoint;

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

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }
            }
        }
    }
    else if(data->mode == DIAG_ORB_POINCARE) {

        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {
            /* Mask dummy markers and those whose time-step was rejected. */
            if( p_f->id[i] > 0 && (p_f->mileage[i] != p_i->mileage[i]) ) {

                real k;
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * data->Npnt + ipoint;

                /* Check and store toroidal crossings. */
                for(int j=0; j < data->ntoroidalplots; j++) {
                    k = diag_orb_check_plane_crossing(p_f->phi[i], p_i->phi[i],
                                                      data->toroidalangles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store poloidal crossings. */
                for(int j=0; j < data->npoloidalplots; j++) {
                    k = diag_orb_check_plane_crossing(p_f->theta[i],
                    		                          p_i->theta[i],
                                                      data->poloidalangles[j]);
                    //printf("Poloidal %10.5f,%10.5f,%10.5d,%10.5d \n",k,data->pncrid[idx],data->npoloidalplots, data->nradialplots);
                    //printf("Poloidal/angles %10.5f,%10.5f,%10.5f,\n,%10.5f,%10.5f\n",data->poloidalangles[j],data->radialdistances[j]
                    //      ,data->rho[imrk * data->Npnt + ipoint],data->r[imrk * data->Npnt + ipoint],data->z[imrk * data->Npnt + ipoint]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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
                        data->pncrid[idx] = j + data->ntoroidalplots;

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                        //printf("Poloidal %10.5f,%10.5f,%10.5f,%10.5f \n",k, p_f->theta[i],
                    	//	p_i->theta[i], data->poloidalangles[j]);;
                    }
                }

                /* Check and store radial crossings. */
                for(int j=0; j < data->nradialplots; j++) {
                    k = diag_orb_check_radial_crossing(p_f->rho[i],p_i->rho[i],
                                                    data->radialdistances[j]);
                    //printf("Radial %10.5f,%10.5f,%10.5f,%10.5f,%10.5d \n",k, p_f->r[i],
                    //		p_i->r[i], data->radialdistances[j],data->nradialplots);
                    //toiseksi viimeinen ei toimi
                    if(k) {
                        real d = k;
			k = 1-d;
                        idx = imrk * data->Npnt + ipoint;
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
                        data->pncrid[idx] = j + data->ntoroidalplots + data->npoloidalplots;

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                        //printf("Radial %10.5f,%10.5f,%10.5f,%10.5f \n",k, p_f->rho[i],
                    	//	p_i->rho[i], data->radialdistances[j]);
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

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + data->writeInterval
                    - p_f->mileage[i];

                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                                    idx = imrk * data->Npnt + ipoint;

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
                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }
            }
        }
    }

    else if(data->mode == DIAG_ORB_POINCARE) {
        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {
            /* Mask dummy markers and thosw whose time-step was rejected. */
            if( p_f->id[i] > 0 && (p_f->mileage[i] != p_i->mileage[i]) ) {

                real k;
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * data->Npnt + ipoint;

                /* Check and store toroidal crossings. */
                for(int j=0; j < data->ntoroidalplots; j++) {
                    k = diag_orb_check_plane_crossing(p_f->phi[i], p_i->phi[i],
                                                      data->toroidalangles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store poloidal crossings. */
                for(int j=0; j < data->npoloidalplots; j++) {
                    k = diag_orb_check_plane_crossing(p_f->theta[i],
                                                      p_i->theta[i],
                                                      data->poloidalangles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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
                        data->pncrid[idx] = j + data->ntoroidalplots;

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                        //printf("Poloidal %10.5f,%10.5f \n",k,data->pncrid[idx]);
                    }
                }

                /* Check and store radial crossings. */
                for(int j=0; j < data->nradialplots; j++) {
                    k = diag_orb_check_radial_crossing(p_f->rho[i],
                                                      p_i->rho[i],
                                                      data->radialdistances[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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
                        		j + data->ntoroidalplots + data->npoloidalplots;

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                        printf("Radial %10.5f,%10.5f \n",k,data->pncrid[idx]);
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

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }

                /* Record marker if enough time has passed from last record, or
                   if marker has met some end condition. */
                real dt = data->mrk_recorded[imrk] + data->writeInterval
                    - p_f->mileage[i];
                if( dt <= 0 || p_f->endcond[i] > 0 ) {
                    idx = imrk * data->Npnt + ipoint;
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

                    ipoint++;
                    if(ipoint == data->Npnt) {
                        ipoint = 0;
                    }
                    data->mrk_pnt[imrk]      = ipoint;
                    data->mrk_recorded[imrk] = p_f->mileage[i];
                }
            }
        }
    }
    else if(data->mode == DIAG_ORB_POINCARE) {
        #pragma omp simd
        for(int i= 0; i < NSIMD; i++) {
            /* Mask dummy markers and thosw whose time-step was rejected. */
            if( p_f->id[i] > 0 && (p_f->mileage[i] != p_i->mileage[i]) ) {

                real k;
                integer imrk   = p_f->index[i];
                integer ipoint = data->mrk_pnt[imrk];
                integer idx    = imrk * data->Npnt + ipoint;

                /* Check and store toroidal crossings. */
                for(int j=0; j < data->ntoroidalplots; j++) {
                    k = diag_orb_check_plane_crossing(p_f->phi[i], p_i->phi[i],
                                                      data->toroidalangles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }

                /* Check and store poloidal crossings. */
                for(int j=0; j < data->npoloidalplots; j++) {
                    k = diag_orb_check_plane_crossing(p_f->theta[i],
                                                      p_i->theta[i],
                                                      data->poloidalangles[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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
                        data->pncrid[idx] = j + data->ntoroidalplots;

                        ipoint++;
                        if(ipoint == data->Npnt) {
                            ipoint = 0;
                        }
                        data->mrk_pnt[imrk]      = ipoint;
                        data->mrk_recorded[imrk] = p_f->mileage[i];
                    }
                }
  
                /* Check and store radial crossings. */
                for(int j=0; j < data->nradialplots; j++) {
                    k = diag_orb_check_radial_crossing(p_f->rho[i],
                                                      p_i->rho[i],
                                                      data->radialdistances[j]);
                    if(k) {
                        real d = 1-k;
                        idx = imrk * data->Npnt + ipoint;
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
			  j + data->ntoroidalplots + data->npoloidalplots;

                        ipoint++;
                        if(ipoint == data->Npnt) {
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
 * Assumption: (ian-fang) < 360 degrees in a single timestep
 *
 * @param fang marker final angle in radians.
 * @param iang marker initial angle in radians.
 * @param ang0 Poincare plane angle in degrees.
 * k = (iang-ang0)/(fang-iang)
 * @return zero if no-crossing, number k, ang0 = k * (fang - iang), otherwise.
 */


real diag_orb_check_plane_crossing(real fang, real iang, real ang0){
  
    real k = 0;
    int fix = iang/CONST_2PI;
    real iaf = fmod(iang,fix);
    real faf = fmod(fang,fix);
    real ang0r = ang0;
    
    if(fix > 0){
        iaf = fmod(iang,CONST_2PI);
	//    printf("iaf %f\n",iaf);
        faf = fmod(fang,CONST_2PI);
	//    printf("faf %f\n",faf);
    }
    
    if(fix < 0){
        iaf = -fix*CONST_2PI + fmod(iang,CONST_2PI);
	//    printf("iaf %f\n",iaf);
        faf = -fix*CONST_2PI + fmod(fang,CONST_2PI);
	//    printf("faf %f\n",faf);
    }
    
    if(iaf < faf){
        if(iaf <= ang0r && faf > ang0r){
	        k = fabs((iaf - ang0r) / (faf - iaf));
        }
        //else if(iaf - CONST_2PI <= ang0 && faf - CONST_2PI > ang0){
        //    k= fabs((iaf - ang0 - CONST_2PI) / (faf - iaf));
        //}
    }
    else{
        if(faf <= ang0r && iaf > ang0r){
            k = fabs((iaf - ang0r) / (faf - iaf));
        }
	//else if(iaf + CONST_2PI <= ang0 && faf + CONST_2PI > ang0){
        //    k = fabs((iaf - ang0 + CONST_2PI) / (faf - iaf));
        //}
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
