#ifndef MCCC_WIENER_H
#define MCCC_WIENER_H

#include "../ascot5.h"

#define MCCC_WIENER_EXCEEDEDCAPACITY 10
#define MCCC_WIENER_NOASSOCIATEDPROCESS 11

/**
*	Struct for storing Wiener processes. Elements of this struct should
*	not be changed (directly) outside mccc package.
*
*	nextslot - integer array where each element shows where the next
*			   wiener process is located (itself if last process). Starts from 0.
*	Nslot    - number of slots, i.e., maximum number of time step reductions
*   Ndim     - Wiener process dimension
*	time     - time instans for different Wiener processes
*	wiener   - Ndim x Nslot array of Wiener process values
*/
typedef struct {
   int* nextslot;
   int Nslot;
   int Ndim;
   real* time;
   real* wiener;
} mccc_wienarr;

mccc_wienarr* mccc_wiener_allocate(int Ndim, int Nslots, real initime);

void mccc_wiener_deallocate(mccc_wienarr* w);

void mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, int* err);

void mccc_wiener_clean(mccc_wienarr* w, real t, int* err);

void mccc_wiener_boxmuller(real* randVar, int Ndim);

void mccc_wiener_error(int err);

#endif
