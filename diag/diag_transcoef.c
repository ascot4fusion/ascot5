/**
 * @file diag_transcoef.c
 * @brief Transport coefficient diagnostics.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../ascot5.h"
#include "../math.h"
#include "../consts.h"
#include "../simulate.h"
#include "../particle.h"
#include "diag_transcoef.h"

#pragma omp declare target
#pragma omp declare simd
real diag_transcoef_check_omp_crossing(real fang, real iang);
#pragma omp declare simd uniform(data)
void diag_transcoef_record(diag_transcoef_data* data, integer index,
                           integer id, real rho, real r, real pitchsign,
                           real t_f, real t_i, real theta_f, real theta_i);
void diag_transcoef_process_and_clean(diag_transcoef_data* data,
                                      integer index, integer id);
#pragma omp end declare target

/**
 * @brief Initializes orbit diagnostics offload data.
 *
 * @param data transport coefficient diagnostics data struct
 * @param offload_data transport coefficient diagnostics offload data struct
 * @param offload_array offload data array
 */
void diag_transcoef_init(diag_transcoef_data* data,
                         diag_transcoef_offload_data* offload_data,
                         real* offload_array) {
    data->id    = &(offload_array[0*offload_data->Nmrk]);
    data->Kcoef = &(offload_array[1*offload_data->Nmrk]);
    data->Dcoef = &(offload_array[2*offload_data->Nmrk]);

    data->interval  = offload_data->interval;
    data->recordrho = offload_data->recordrho;
    data->Navg      = offload_data->Navg;

    data->datapoints = malloc(
        offload_data->Nmrk*sizeof(diag_transcoef_link*) );
    for(int i = 0; i < offload_data->Nmrk; i++) {
        data->id[i] = -1;
        data->datapoints[i] = NULL;
    }
}

void diag_transcoef_free(diag_transcoef_data* data) {
    free(data->datapoints);

}

/**
 * @brief Collect transport diagnostics for fo simulation.
 *
 * @param data transport diagnostics data struct
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void diag_transcoef_update_fo(diag_transcoef_data* data,
                              particle_simd_fo* p_f, particle_simd_fo* p_i) {
    #pragma omp simd
    for(int i=0; i < NSIMD; i++) {
        real p[3] = {p_f->p_r[i], p_f->p_phi[i], p_f->p_z[i]};
        real B[3] = {p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]};
        real pitchsign = 1 - 2*(math_dot(p, B) < 0);
        diag_transcoef_record(
            data, p_f->index[i], p_f->id[i], p_f->rho[i], p_f->r[i], pitchsign,
            p_f->mileage[i], p_i->mileage[i], p_f->theta[i],  p_i->theta[i]);
    }

    /* If marker simulation was ended, process and clean the data */
    for(int i=0; i < NSIMD; i++) {

        /* Mask dummy markers and those which are running */
        if( p_f->id[i] < 1 || p_f->running[i] > 0 ) {
            continue;
        }

        diag_transcoef_process_and_clean(data, p_f->index[i], p_f->id[i]);

    }
}

/**
 * @brief Collect transport diagnostics for gc simulation.
 *
 * @param data transport diagnostics data struct
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void diag_transcoef_update_gc(diag_transcoef_data* data,
                              particle_simd_gc* p_f, particle_simd_gc* p_i) {
    #pragma omp simd
    for(int i=0; i < NSIMD; i++) {
        real pitchsign = 1 - 2*(p_f->ppar[i] < 0);
        diag_transcoef_record(
            data, p_f->index[i], p_f->id[i], p_f->rho[i], p_f->r[i], pitchsign,
            p_f->mileage[i], p_i->mileage[i], p_f->theta[i],  p_i->theta[i]);
    }


    /* If marker simulation was ended, process and clean the data */
    for(int i=0; i < NSIMD; i++) {

        /* Mask dummy markers and those which are running */
        if( p_f->id[i] < 1 || p_f->running[i] > 0 ) {
            continue;
        }

        diag_transcoef_process_and_clean(data, p_f->index[i], p_f->id[i]);
    }
}

/**
 * @brief Collect transport diagnostics for ml simulation.
 *
 * @param data transport diagnostics data struct
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void diag_transcoef_update_ml(diag_transcoef_data* data,
                              particle_simd_ml* p_f, particle_simd_ml* p_i) {
    #pragma omp simd
    for(int i=0; i < NSIMD; i++) {
        real pitchsign = 1 - 2*(p_f->pitch[i] < 0);
        diag_transcoef_record(
            data, p_f->index[i], p_f->id[i], p_f->rho[i], p_f->r[i], pitchsign,
            p_f->mileage[i], p_i->mileage[i], p_f->theta[i],  p_i->theta[i]);
    }


    /* If marker simulation was ended, process and clean the data */
    for(int i=0; i < NSIMD; i++) {
        /* Mask dummy markers and those which are running */
        if( p_f->id[i] < 1 || p_f->running[i] > 0 ) {
            continue;
        }

        diag_transcoef_process_and_clean(data, p_f->index[i], p_f->id[i]);
    }
}


/**
 * @brief Check if criteria for recording is met for a single marker and record.
 *
 * @param data pointer to transport coefficient data.
 * @param index marker index in the marker queue.
 * @param id marker id.
 * @param rho.
 * @param r.
 * @param pitchsign.
 * @param t_f time (mileage) at the end of the time step.
 * @param t_i time (mileage) at the beginning of the time step.
 * @param theta_f poloidal angle at the end of the time step.
 * @param theta_i poloidal angle at the beginning of the time step.
 */
void diag_transcoef_record(diag_transcoef_data* data, integer index,
                           integer id, real rho, real r, real pitchsign,
                           real t_f, real t_i, real theta_f, real theta_i) {
    /* Mask dummy markers */
    if( id > 0 ) {

        /* Check whether marker position should be recorded: *
         * - Time step was accepted t_f > t_i
         * - Enough time has passed since last record t_f - ti > interval OR
             no data exists yet.
         * - Marker has crossed OMP during the current time step.
         */
        real record = 0.0;
        if( t_f > t_i ) {
            if( data->datapoints[index] == NULL ) {
                record = diag_transcoef_check_omp_crossing(theta_f, theta_i);
            }
            else if( t_f - data->datapoints[index]->time > data->interval ) {
                record = diag_transcoef_check_omp_crossing(theta_f, theta_i);
            }
        }

        /* Record */
        if( record > 0) {
            diag_transcoef_link* newlink =
                malloc(sizeof(diag_transcoef_link));
            newlink->time      = t_f;
            newlink->pitchsign = pitchsign;

            if(data->recordrho) {
                newlink->rho   = rho;
            }
            else {
                newlink->rho   = r;
            }

            /* Store the link or make a new chain if the marker is new */
            if( data->datapoints[index] == NULL ) {
                newlink->prevlink   = NULL;
                data->datapoints[index] = newlink;
            }
            else {
                newlink->prevlink   = data->datapoints[index];
                data->datapoints[index] = newlink;
            }
        }
    }
}


/**
 * @brief Process recorded data to transport coefficients and clean.
 *
 * This function is called when marker simulation has ended.
 *
 * @param data pointer to transport coefficient data.
 * @param index marker index in the marker queue.
 * @param id marker id.
 * @param rho.
 * @param r.
 * @param pitchsign.
 * @param t_i time (mileage) at the beginning of the time step.
 * @param t_f time (mileage) at the end of the time step.
 * @param theta_i poloidal angle at the beginning of the time step.
 * @param theta_f poloidal angle at the end of the time step.
 */
void diag_transcoef_process_and_clean(diag_transcoef_data* data,
                                      integer index, integer id) {

    /* Count number of positive and negative crossings */
    int positive = 0, negative = 0;
    diag_transcoef_link* link = data->datapoints[index];
    while(link != NULL) {
        if(link->pitchsign < 0) {
            negative++;
        }
        else {
            positive++;
        }
        link = link->prevlink;
    }

    /* Which ever there are more are stored */
    int datasize, sign;
    if(positive >= negative) {
        datasize = positive;
        sign     = 1;
    }
    else {
        datasize = negative;
        sign     = -1;
    }

    /* If there are enough datapoints, process them to K and D */
    if(datasize > data->Navg) {
        /* How many points we have after averaging data */
        int navgpnt = ceil(datasize/data->Navg);
        real* rho  = malloc(navgpnt*sizeof(real));
        real* time = malloc(navgpnt*sizeof(real));

        /* The datasize is not necessarily multiple of navg. We take the
           "extra" points from last points and average them separately   */
        link = data->datapoints[index];
        while(link->pitchsign * sign < 0) {
            link = link->prevlink;
        }
        rho[navgpnt-1]  = link->rho;
        time[navgpnt-1] = link->time;

        int nextrapnts = datasize - navgpnt*data->Navg;
        if(nextrapnts == 0) {
            nextrapnts = data->Navg;
        }
        for(int k = 1; k < nextrapnts; k++) {
            link = link->prevlink;
            while(link->pitchsign * sign < 0) {
                link = link->prevlink;
            }
            rho[navgpnt-1]  += link->rho;
            time[navgpnt-1] += link->time;
        }
        rho[navgpnt-1]  /= nextrapnts;
        time[navgpnt-1] /= nextrapnts;

        /* Take the full averages */
        for(int j = navgpnt-2; j > -1; j--) {
            rho[j]  = 0;
            time[j] = 0;
            for(int k = 0; k < data->Navg; k++) {
                link = link->prevlink;
                while(link->pitchsign * sign < 0) {
                    link = link->prevlink;
                }
                rho[j]  += link->rho;
                time[j] += link->time;
            }
            rho[j]  /= data->Navg;
            time[j] /= data->Navg;
        }

        /* Evaluate coefficients */
        real K = 0;
        for(int j = 0; j < navgpnt-1; j++) {
            K += ( rho[j+1] - rho[j] ) / ( time[j+1] - time[j] );
        }
        K = K / navgpnt;

        real D = 0;
        for(int j = 0; j < navgpnt-1; j++) {
            real a = rho[j+1] - rho[j] - K * ( time[j+1] - time[j] );
            D += 0.5*a*a / ( time[j+1] - time[j] );
        }
        D = D / navgpnt;

        /* Store and free resources */
        free(rho);
        free(time);
        data->id[index]    = (real)id;
        data->Kcoef[index] = K;
        data->Dcoef[index] = D;
    }

    /* Clear temporary storage */
    diag_transcoef_link* temp;
    link = data->datapoints[index];
    while(link != NULL) {
        temp = link->prevlink;
        free(link);
        link = temp;
    }
    data->datapoints[index] = NULL;

}


/**
 * @brief Check if marker has crossed omp.
 *
 * @param fang marker final poloidal angle in radians.
 * @param iang marker initial poloidal angle in radians.
 *
 * @return zero if no-crossing, number k, k = (iang - fang), otherwise.
 */
real diag_transcoef_check_omp_crossing(real fang, real iang){

    real k = 0.0;
    if( floor( fang/CONST_2PI ) != floor( iang/CONST_2PI ) ) {

        real a = fmod(iang, CONST_2PI);
        if(a < 0){
            a = CONST_2PI + a;
        }

        a = fabs(a);
        if(a > CONST_PI) {
            a = CONST_2PI - a;
        }
        k = fabs(a / (fang - iang));
    }

    return k;
}
