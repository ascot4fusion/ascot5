/** @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
*
* A module for generating Wiener processes. When adaptive time 
* step is used (and steps are rejected), Wiener processes are
* generated using the so-called Brownian bridge. This module
* contains associated helper routines as well.
*/
#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "math.h"
#include "../ascot5.h"
#include "../consts.h"
#include "mccc_wiener.h"

int MCCC_EMPTY = -999;

/**	Allocates an array that will be used to store generated Wiener processes until
*	they are no longer needed. Deallocation is done with subroutine wiener_deallocate.
*	Array has format (Ndim+2) X Nslots where the first row indicates the position of 
*	the next process in time or its own position if no such process exists. The first
*	process is always located at the firs column. The second row stores the time 
*	instants processes correspond. The rest of the rows stores values of the Ndim
*	dimensional Wiener processes. Nslots value should not be overly large.
*
*	integer Ndim         -- Number of Wiener process dimensions
*	integer Nslots       -- Maximum number of time instances for which processes are stored.
*	real*8  initime      -- Time (>0) at the beginning of the simulation [s]
*	real*8  wienarr(:,:) -- array that stores the Wiener processes
*/
mccc_wienarr* mccc_wiener_allocate(int Ndim, int Nslots, real initime){

    int i;

    mccc_wienarr* w = malloc(sizeof(mccc_wienarr));
    w->Nslot = Nslots;
    w->Ndim = Ndim;
    
    w->wiener = malloc(sizeof(real)*Ndim*Nslots);
    w->time = malloc(sizeof(real)*Nslots);
    w->nextslot = malloc(sizeof(int)*Nslots);

    // Initialize position instances indicating all slots are empty
    for(i = 0; i < Nslots; i = i +1){
	w->time[i] = MCCC_EMPTY;
	w->nextslot[i] = MCCC_EMPTY;
    }

    // W(t_0) = 0 by the definition of the Wiener process. Here we initialize it.
    w->nextslot[0] = 0;
    w->time[0] = initime;
    
    for(i = 0; i < Ndim; i = i +1){
       w->wiener[i] = 0.0;
    }

    return w;
}


/** Deallocates the array storing the Wiener processes.
*
*	w struct that stores the Wiener processes
*/
void mccc_wiener_deallocate(mccc_wienarr* w){
    free(w->nextslot);
    free(w->time);
    free(w->wiener);
    free(w);
}


/** Generates a new Wiener process at a given time instant taking the Brownian bridge
* 	into account if needed.
*
*	@param w      array that stores the Wiener processes
*	@param t      time for which the new process will be generated
*	@param windex index of the generated Wiener process in the Wiener array
*	@param err    error flag, negative indicates something went wrong
*/
void mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, int* err){
    
    int idx, eidx, i; // Helper variables
    int im, ip; // Indexes of the Wiener processes for which tm < t < tp

    int Nslots = w->Nslot;
    int Ndim = w->Ndim; 
    *windex = -1;
    *err=0;

    ip = -1; // There isn't necessarily a Wiener process tp > t

    // Find im and ip
    idx = 0;
    for(i = 0; i < Nslots; i=i+1){
	if(w->time[idx] < t) {
	    im = idx;
	}
	else if(w->time[idx] == t){
	    // It seems that the process we are looking for already exists!
	    *windex = idx;
	    return;
	}
	else if(w->time[idx] > t){
	    ip = idx;
	    break;
	}
	
	if(w->nextslot[idx] == idx){
	    // Reached last process
	    break;
	}
	
	idx = w->nextslot[idx];
    }
    
    // Find an empty slot for the next process
    eidx = 0;
    for( i = 0; i < Nslots; i = i+1 ){
	if( w->nextslot[i] == MCCC_EMPTY){
	    eidx = i;
	    break;
	}
    }
    if(eidx == 0){
	// It seems that we have exceeded capacity of the Wiener array
	// Produce an error.
	*err = MCCC_WIENER_EXCEEDEDCAPACITY;
	return;
    }

    
    // The eidx entry in the wiener array is always empty. We use that for temporary storage
    // to spare one allocation/deallocation cycle.
    mccc_wiener_boxmuller( &(w->wiener[(Ndim)*eidx]), Ndim);
    if(ip == -1){
	// There are no Wiener processes existing for tp > t.
	// The generated Wiener process then has a mean W(tm) and variance t-tm.

	w->nextslot[eidx] = eidx;
	w->time[eidx] = t;
	for(i=0; i < Ndim; i=i+1){
	    w->wiener[i + eidx*w->Ndim] = w->wiener[i + im*w->Ndim] + sqrt(t-w->time[im])*w->wiener[i + eidx*w->Ndim];
	}
	*windex = eidx;
	w->nextslot[im] = eidx;
	}
    else{
	// A Wiener process for tp > t exist. Generate a new process using the rules
	// set by the Brownian bridge. The rules are:
	//
	// mean = W(tm) + ( W(ip)-W(im) )*(t-tm)/(tp-tm)
	// variance = (t-tm)*(tp-t)/(tp-tm)
	w->time[eidx] = t;
	for(i=0;i < Ndim; i = i+1){
	    w->wiener[i + eidx*w->Ndim] = w->wiener[i + im*w->Ndim] + ( w->wiener[i + ip*w->Ndim] - w->wiener[i + im*w->Ndim] )
		*( t-w->time[im] )/( w->time[ip]-w->time[im] )
		+ sqrt( ( t-w->time[im] )*( w->time[ip]-t )/( w->time[ip]-w->time[im] ) )
		*w->wiener[i + eidx*w->Ndim];
	}
	// Sort new wiener process to its correct place
	w->nextslot[eidx] = ip;
	w->nextslot[im] = eidx;
	*windex = eidx;
    }
}

/**	Removes Wiener processes from the array that are no longer required. Processes W(t')  
*	are redundant if t' <  t, where t is the current simulation time. Note that W(t) 
*	should exist before W(t') are removed. This routine should be called each time when 
*	simulation time is advanced.
*	
*	@param w   array that stores the Wiener processes
*	@param t   time for which the new process will be generated
*	@param err error flag, negative indicates something went wrong
*/
 void mccc_wiener_clean(mccc_wienarr* w, real t, int* err){

     int idx, nidx, i, nextidx;

     *err = 0;
        
     // Remove processes W(t_i) until ti = t
     idx = 0;
     real ti = w->time[idx];
     while(ti < t){
	 nextidx = w->nextslot[idx];
	 if(idx == nextidx){
	     *err = MCCC_WIENER_NOASSOCIATEDPROCESS;
	     return;
	 }
	 
	 w->nextslot[idx] = MCCC_EMPTY;
	 w->time[idx] = MCCC_EMPTY;
	 
	 idx = nextidx;
	 ti = w->time[idx];
     }

     if(idx==0){
	 return;
     }
     
     // Move W(t) process as the first one
     w->nextslot[0] = w->nextslot[idx];
     
     w->time[0] = w->time[idx];
     w->time[idx] = MCCC_EMPTY;
     for(i = 0; i < w->Ndim; i=i+1){
	 w->wiener[i] = w->wiener[idx*w->Ndim+i];
     }
     
     // Check if the process is also the last one
     if( w->nextslot[idx] == idx ){
	 w->nextslot[0] = 0;
     }
     w->nextslot[idx] = MCCC_EMPTY;

}
  
  

/**	Generates standard normally distributed random numbers using the geometric Box-Muller
* 	method.
*
* 	@param randVar array to be populated with random numbers
*	@param Ndim dimension of the array
*/
void mccc_wiener_boxmuller(real* randVar, int Ndim){
    
    real x1, x2, w, s;
    int isOdd; // Indicates if even number of random numbers are requested
    int i; // Helper variables
    
    isOdd = (Ndim+1) % 2;
#ifdef MCCC_WIENER_USE_GBM
    // The geometric form
    for( i = 0; i < Ndim; i=i+2){
	w = 2.0;
	while( w >= 1.0 ){
	    x1 = 2*((real)drand48())-1;
	    x2 = 2*((real)drand48())-1;
	    w = x1*x1 + x2*x2;
	}
	
	w = sqrt( (-2 * log( w ) ) / w );
	randVar[i] = x1 * w;
	if((i < Ndim-2) || (isOdd > 0)) {
	    randVar[i+1] = x2 * w;
	}
    }
#else
    // The common form which could be faster on some platforms
    for( i = 0; i < Ndim; i=i+2){
	x1 = ((real)drand48());
	x2 = (drand48());
	w = sqrt(-2*log(x1));
	s = cos(CONST_2PI*x2);
	randVar[i] = w*s;
	if((i < Ndim-2) || (isOdd > 0) ){
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

