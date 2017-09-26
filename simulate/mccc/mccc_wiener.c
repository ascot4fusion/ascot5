/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file mccc_wiener.c
 * @brief A module for handling Wiener processes
 *
 * A module for handling Wiener processes. When adaptive time 
 * step is used (and steps are rejected), Wiener processes are
 * generated using the so-called Brownian bridge. This module
 * contains associated helper routines.
 */
#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../math.h"
#include "../../ascot5.h"
#include "../../consts.h"
#include "mccc_wiener.h"

const int MCCC_EMPTY = -999;

/**
 * @brief Initializes a struct that will be used to store generated Wiener processes 
 *
 * @param w Wiener struct to be initialized
 * @param initime time instance corresponding to initial Wiener process (which has value W(t) = 0)
 */
void mccc_wiener_initialize(mccc_wienarr* w, real initime){

    int i;

    /* Initialize position instances indicating all slots are empty */
    for(i = 0; i < MCCC_NSLOTS; i = i +1){
	w->time[i] = MCCC_EMPTY;
	w->nextslot[i] = MCCC_EMPTY;
    }

    /* W(t_0) = 0 by the definition of the Wiener process. Here we initialize it. */
    w->nextslot[0] = 0;
    w->time[0] = initime;
    
    for(i = 0; i < MCCC_NDIM; i = i +1){
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
 * @param err error flag, negative indicates something went wrong
 */
void mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, real* rand5, int* err){
    
    int eidx, i; /* Helper variables */
    int im, ip; /* Indexes of the Wiener processes for which tm < t < tp */

    *windex = -1;
    *err=0;

    ip = -1; /* There isn't necessarily a Wiener process tp > t */

    /* Find im and ip */
    int idx = 0;
    im = 0;
    for(i = 0; i < MCCC_NSLOTS; i=i+1){
	if(w->nextslot[idx] == idx){
	    /* Reached last process */
	    i = MCCC_NSLOTS;
	}
        
	if(w->time[idx] == t) {
	    /* Process already exists */
	    *windex = idx;
	    i = MCCC_NSLOTS;
	}
	else {
	    if(w->time[idx] < t ) {
		/* Process i for which t_i < t */
		im = idx;
	    }
	    if(w->time[w->nextslot[idx]] > t) {
		/* Process i for which t_i > t */
		ip = w->nextslot[idx];
		i = MCCC_NSLOTS; // Exit loop
	    }
	}
	    
	idx = w->nextslot[idx];
	
    }
    
    if(*windex == -1) {
	/* Find an empty slot for the next process */
	eidx = 0;
	for( i = 0; i < MCCC_NSLOTS; i = i+1 ){
	    if( w->nextslot[i] == MCCC_EMPTY){
		eidx = i;
		i = MCCC_NSLOTS;
	    }
	}
	if(eidx == 0){
	    /* It seems that we have exceeded capacity of the Wiener array
	     * Produce an error. */
	    *err = MCCC_WIENER_EXCEEDEDCAPACITY;
	}

	if(!(*err)) {

	    if(ip == -1){
		/* There are no Wiener processes existing for tp > t.
		 * The generated Wiener process then has a mean W(tm) and variance t-tm. */

		w->nextslot[eidx] = eidx;
		w->time[eidx] = t;
		for(i=0; i < MCCC_NDIM; i=i+1){
		    w->wiener[i + eidx*MCCC_NDIM] = w->wiener[i + im*MCCC_NDIM] + sqrt(t-w->time[im])*rand5[i];
		}
		*windex = eidx;
		w->nextslot[im] = eidx;
	    }
	    else{
		/* A Wiener process for tp > t exist. Generate a new process using the rules
		 * set by the Brownian bridge. The rules are:
		 *
		 * mean = W(tm) + ( W(ip)-W(im) )*(t-tm)/(tp-tm)
		 * variance = (t-tm)*(tp-t)/(tp-tm) */
		w->time[eidx] = t;
		for(i=0;i < MCCC_NDIM; i = i+1){
		    w->wiener[i + eidx*MCCC_NDIM] = w->wiener[i + im*MCCC_NDIM] + ( w->wiener[i + ip*MCCC_NDIM] - w->wiener[i + im*MCCC_NDIM] )
			*( t-w->time[im] )/( w->time[ip]-w->time[im] )
			+ sqrt( ( t-w->time[im] )*( w->time[ip]-t )/( w->time[ip]-w->time[im] ) )
			*rand5[i];
		}
		/* Sort new wiener process to its correct place */
		w->nextslot[eidx] = ip;
		w->nextslot[im] = eidx;
		*windex = eidx;
	    }
	}
    }
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
 * @param err error flag, negative indicates something went wrong
 */
void mccc_wiener_clean(mccc_wienarr* w, real t, int* err){

    int idx, nidx, i, nextidx;

    *err = 0;
        
    /* Remove processes W(t_i) until ti = t */
    idx = 0;
    real ti = w->time[idx];
    while(ti < t){
	nextidx = w->nextslot[idx];
	if(idx == nextidx){
	    *err = MCCC_WIENER_NOASSOCIATEDPROCESS;
	    t = ti; // Breaks the loop
	}
	else {
	    w->nextslot[idx] = MCCC_EMPTY;
	    w->time[idx] = MCCC_EMPTY;
	 
	    idx = nextidx;
	    ti = w->time[idx];
	}
    }

    if(idx!=0 && !(*err)){
	/* Move W(t) process as the first one */
	w->nextslot[0] = w->nextslot[idx];
     
	w->time[0] = w->time[idx];
	w->time[idx] = MCCC_EMPTY;
	for(i = 0; i < MCCC_NDIM; i=i+1){
	    w->wiener[i] = w->wiener[idx*MCCC_NDIM+i];
	}
     
	/* Check if the process is also the last one */
	if( w->nextslot[idx] == idx ){
	    w->nextslot[0] = 0;
	}
	w->nextslot[idx] = MCCC_EMPTY;
    }

}
  
  
/**
 * @brief Generates standard normally distributed random numbers 
 * 
 * Random numbers are created using the Box-Muller method.
 *
 * Compiler flag  A5_CCOL_USE_GEOBM in ascot5.h determines whether geometric
 * or common form is used.
 *
 * @param randVar pointer to array to be populated with random numbers
 * @param Ndim dimension of the array
 *
 * @todo Move to math.h
 * @todo Implement MCCC_WIENER_USE_GBM flag
 */
void mccc_wiener_boxmuller(real* randVar, int Ndim){
    
    real x1, x2, w, s;
    int isEven; /* Indicates if even number of random numbers are requested */
    int i; /* Helper variables */
    
    isEven = (Ndim+1) % 2;
#if A5_CCOL_USE_GEOBM == 1
    /* The geometric form */
    for( i = 0; i < Ndim; i=i+2){
	w = 2.0;
	while( w >= 1.0 ){
	    x1 = 2*((real)drand48())-1;
	    x2 = 2*((real)drand48())-1;
	    w = x1*x1 + x2*x2;
	}
	
	w = sqrt( (-2 * log( w ) ) / w );
	randVar[i] = x1 * w;
	if((i < Ndim-2) || (isEven > 0)) {
	    randVar[i+1] = x2 * w;
	}
    }
#else
    /* The common form */
    for( i = 0; i < Ndim; i=i+2){
	x1 = ((real)drand48());
	x2 = (drand48());
	w = sqrt(-2*log(x1));
	s = cos(CONST_2PI*x2);
	randVar[i] = w*s;
	if((i < Ndim-2) || (isEven > 0) ){
	    if(x2 < 0.5){
		randVar[i+1] = w*sqrt(1-s*s);
	    }
	    else{
		randVar[i+1] = -w*sqrt(1-s*s);
	    }
	}
    }
#endif

}

