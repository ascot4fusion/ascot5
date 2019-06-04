/**
 * @file diag_transcoef.c
 * @brief Transport coefficient diagnostics.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../simulate.h"
#include "../particle.h"
#include "diag_transcoef.h"

#pragma omp declare target
#pragma omp declare simd
real diag_transcoef_check_omp_crossing(real fang, real iang);
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

    data->interval = offload_data->interval;
    data->Navg     = offload_data->Navg;

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

void diag_transcoef_update_gc(diag_transcoef_data* data,
                              particle_simd_gc* p_f, particle_simd_gc* p_i) {
    #pragma omp simd
    for(int i=0; i < NSIMD; i++) {
        /* Mask dummy markers */
        if( p_f->id[i] > 0 ) {

            /* Check whether marker positin should be recorded */
            real record = 0.0;
            if( data->interval >= 0 &&
                (p_f->time[i] - p_i->time[i]) > data->interval ) {
                record = (p_f->time[i] - p_i->time[i]) / data->interval;
            }
            if( data->interval < 0 && p_f->time[i] > p_i->time[i] ) {
                record = diag_transcoef_check_omp_crossing(
                    p_f->theta[i], p_i->theta[i]);
            }

            /* Record */
            if( record > 0) {
                diag_transcoef_link* newlink =
                    malloc(sizeof(diag_transcoef_link));
                newlink->rho       = p_f->rho[i];
                newlink->time      = p_f->time[i];
                newlink->pitchsign = 1 - 2*(p_f->vpar[i] < 0);

                /* Store the link or make a new chain if the marker is new */
                if( data->datapoints[p_f->index[i]] == NULL ) {
                    newlink->prevlink   = NULL;
                    data->datapoints[p_f->index[i]] = newlink;
                }
                else {
                    newlink->prevlink   = data->datapoints[p_f->index[i]];
                    data->datapoints[p_f->index[i]] = newlink;
                }
            }
        }
    }


    /* If marker simulation was ended, process and clean the data */
    for(int i=0; i < NSIMD; i++) {
        /* Mask dummy markers and those which are running */
        if( p_f->id[i] < 1 || p_f->running[i] > 0 ) {
            continue;
        }

        /* Count number of positive and negative crossings */
        int positive = 0, negative = 0;
        diag_transcoef_link* link = data->datapoints[p_f->index[i]];
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
            link = data->datapoints[p_f->index[i]];
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
            data->id[p_f->index[i]]    = (real)p_f->id[i];
            data->Kcoef[p_f->index[i]] = K;
            data->Dcoef[p_f->index[i]] = D;
        }

        /* Clear temporary storage */
        diag_transcoef_link* temp;
        link = data->datapoints[p_f->index[i]];
        while(link != NULL) {
            temp = link->prevlink;
            free(link);
            link = temp;
        }
        data->datapoints[p_f->index[i]] = NULL;
    }
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

    real ang0 = 0;
    real k = 0.0;
    if( floor( (fang + ang0)/CONST_2PI ) != floor( (iang + ang0)/CONST_2PI ) ) {

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
