/**
 * @file mccc_wiener.c
 * @brief A module for handling Wiener processes
 *
 * A module for handling Wiener processes. When adaptive time step is used (and
 * steps are rejected), Wiener processes are generated using the so-called
 * Brownian bridge. This module contains associated helper routines.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils/mathlib.h"
#include "defines.h"
#include "consts.h"
#include "coulomb_collisions.h"

/** Indicates an empty slot in wiener array */
#ifdef _OPENACC
const int MCCC_EMPTY = -1;
#pragma acc declare copyin(MCCC_EMPTY)
#elif defined(_OPENMP)
DECLARE_TARGET
const int MCCC_EMPTY = -1;
DECLARE_TARGET_END
#else
const int MCCC_EMPTY = -1;
#endif

/**
 * @brief Initializes a struct that stores generated Wiener processes
 *
 * @param w Wiener struct to be initialized
 * @param initime time when a Wiener process begins
 */
void mccc_wiener_initialize(mccc_wienarr* w, real initime){

    /* Initialize position instances indicating all slots are empty */
    for(int i=0; i<MCCC_NSLOTS; i++){
        w->nextslot[i] = MCCC_EMPTY;
    }

    /* W(t_0) = 0 by the definition of the Wiener process */
    w->nextslot[0] = 0;
    w->time[0]     = initime;
    for(int i = 0; i < MCCC_NDIM; i++){
        w->wiener[i] = 0.0;
    }
}

/**
 * @brief Generates a new Wiener process at a given time instant
 *
 * Generates a new Wiener process. The generated process is drawn from
 * normal distribution unless there exists a Wiener process at future
 * time-instance, in which case the process is created using the Brownian
 * bridge.
 *
 * @param w array that stores the Wiener processes
 * @param t time for which the new process will be generated
 * @param windex index of the generated Wiener process in the Wiener array
 * @param rand5 array of 5 normal distributed random numbers
 *
 * @return zero if generation succeeded
 */
err_t mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, real* rand5){
    err_t err = 0;
    int eidx; /* Helper variables */
    int im, ip; /* Indices of the Wiener processes for which tm < t < tp */

    windex[0] = -1;
    im = 0;
    ip = -1; /* There isn't necessarily a Wiener process tp > t */

    /* Find im and ip */
    int idx = 0, breakloop = 0;
    for(int i=0; i<MCCC_NSLOTS; i++){
        if(!breakloop) {
            if(w->nextslot[idx] == idx){
                /* Reached last process. Break loop */
                breakloop = 1;
            }
            if(w->time[idx] == t) {
                /* Process already exists */
                windex[0] = idx;
                breakloop = 1;
            }
            else {
                if(w->time[idx] < t ) {
                    /* Process i for which t_i < t */
                    im = idx;
                }
                if(w->time[w->nextslot[idx]] > t) {
                    /* Process i for which t_i > t */
                    ip = w->nextslot[idx];
                    breakloop = 1;
                }
            }
            idx = w->nextslot[idx];
        }
    }

    if(windex[0] == -1) {
        /* Find an empty slot for the next process */
        eidx = 0;
        for(int i=0; i<MCCC_NSLOTS; i++){
            if( w->nextslot[i] == MCCC_EMPTY){
                eidx = i;
                i = MCCC_NSLOTS;
            }
        }
        if(eidx == 0){
            /* It seems that we have exceeded capacity of the Wiener array
             * Produce an error. */
            err = error_raise(ERR_WIENER_ARRAY, __LINE__, EF_MCCC_WIENER);
        }

        if(!err) {
            if(ip == -1){
                /* There are no Wiener processes existing for tp > t.     *
                 * The generated Wiener process then has a mean W(tm) and *
                 * variance t-tm.                                         */

                w->nextslot[eidx] = eidx;
                w->time[eidx]     = t;
                for(int i=0; i<MCCC_NDIM; i++){
                    w->wiener[i + eidx*MCCC_NDIM] = w->wiener[i + im*MCCC_NDIM]
                        + sqrt(t-w->time[im])*rand5[i];
                }
                windex[0]       = eidx;
                w->nextslot[im] = eidx;
            }
            else{
                /* A Wiener process for tp > t exist. Generate a new process
                 * using the rules set by the Brownian bridge. The rules are:
                 *
                 * mean = W(tm) + ( W(ip)-W(im) )*(t-tm)/(tp-tm)
                 * variance = (t-tm)*(tp-t)/(tp-tm) */
                w->time[eidx] = t;
            for(int i=0; i<MCCC_NDIM; i++){
                w->wiener[i+eidx*MCCC_NDIM] =
                      w->wiener[i+im*MCCC_NDIM]
                    + (   w->wiener[i + ip*MCCC_NDIM]
                        - w->wiener[i + im*MCCC_NDIM] )
                    * ( t - w->time[im] ) / ( w->time[ip] - w->time[im] )
                    + sqrt( ( t - w->time[im] ) * ( w->time[ip] - t )
                            / ( w->time[ip] - w->time[im] ) ) * rand5[i];
            }

            /* Sort new wiener process to its correct place */
            w->nextslot[eidx] = ip;
            w->nextslot[im]   = eidx;
            windex[0]         = eidx;
            }
        }
    }

    return err;
}

/**
 * @brief Removes Wiener processes from the array that are no longer required.
 *
 * Processes W(t') are redundant if t' <  t, where t is the current simulation
 * time. Note that W(t) should exist before W(t') are removed. This routine
 * should be called each time when simulation time is advanced.
 *
 * @param w array that stores the Wiener processes
 * @param t time for which the new process will be generated
 *
 * @return zero if cleaning succeeded
 */
err_t mccc_wiener_clean(mccc_wienarr* w, real t){
    err_t err = 0;
    int idx, nextidx;

    /* Remove processes W(t_i) until ti = t */
    idx = 0;
    real ti = w->time[idx];
    while(ti < t){
        nextidx = w->nextslot[idx];
        if(idx == nextidx){
            err = error_raise(ERR_WIENER_ARRAY, __LINE__, EF_MCCC_WIENER);
            t = ti; // Breaks the loop
        }
        else {
            w->nextslot[idx] = MCCC_EMPTY;
            idx = nextidx;
            ti  = w->time[idx];
        }
    }

    if(idx!=0 && !err){
        /* Move W(t) process as the first one */
        w->nextslot[0] = w->nextslot[idx];

        w->time[0]   = w->time[idx];
        for(int i=0; i<MCCC_NDIM; i++){
            w->wiener[i] = w->wiener[idx*MCCC_NDIM+i];
        }

        /* Check if the process is also the last one */
        if( w->nextslot[idx] == idx ){
            w->nextslot[0] = 0;
        }
        w->nextslot[idx] = MCCC_EMPTY;
    }

    return err;
}
