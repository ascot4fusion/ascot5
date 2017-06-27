/**
 * @file mccc_wiener.h
 * @brief header file for mccc_wiener.c
 */
#ifndef MCCC_WIENER_H
#define MCCC_WIENER_H

#include "../../ascot5.h"

#define MCCC_WIENER_EXCEEDEDCAPACITY 10
#define MCCC_WIENER_NOASSOCIATEDPROCESS 11

/**
 * @brief Struct for storing Wiener processes. Elements of this struct should
 *	  not be changed (directly) outside mccc package.
 */
typedef struct {
    int* nextslot; /** Integer array where each element shows where the next
		    *  wiener process is located (itself if last process). 
		    *  Starts from 0.*/
    int Nslot; /** Number of slots, i.e., maximum number of time step reductions */
    int Ndim; /** Wiener process dimension (5 for guiding centers) */
    real* time; /** Time instances for different Wiener processes */
    real* wiener; /** Ndim x Nslot 1D array of Wiener process values */
} mccc_wienarr;

#pragma omp declare target
#pragma omp declare simd
mccc_wienarr* mccc_wiener_allocate(int Ndim, int Nslots, real initime);

#pragma omp declare simd
void mccc_wiener_deallocate(mccc_wienarr* w);

#pragma omp declare simd
void mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, int* err);

#pragma omp declare simd
void mccc_wiener_clean(mccc_wienarr* w, real t, int* err);

#pragma omp declare simd
void mccc_wiener_boxmuller(real* randVar, int Ndim);

#pragma omp declare simd
void mccc_wiener_error(int err);

#pragma omp end declare target

#endif
