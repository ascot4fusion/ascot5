/**
 * @file mccc_wiener.h
 * @brief header file for mccc_wiener.c
 */
#ifndef MCCC_WIENER_H
#define MCCC_WIENER_H

#include "../../ascot5.h"
#include "../../error.h"

#define MCCC_NDIM 5 /* Wiener process dimension. NDIM=5 because 
			    only guiding centers are simulated with 
			    adaptive time step */
#define MCCC_NSLOTS WIENERSLOTS /* Maximum time step reductions */

/**
 * @brief Struct for storing Wiener processes. Elements of this struct should
 *	  not be changed (directly) outside mccc package.
 */
typedef struct {
    int nextslot[MCCC_NSLOTS]; /** Integer array where each element shows where the next
		    *  wiener process is located (itself if last process). 
		    *  Starts from 0.*/
    real time[MCCC_NSLOTS]; /** Time instances for different Wiener processes */
    real wiener[MCCC_NDIM*MCCC_NSLOTS]; /** Ndim x Nslot 1D array of Wiener process values */
} mccc_wienarr;

#pragma omp declare target
#pragma omp declare simd
void mccc_wiener_initialize(mccc_wienarr* w, real initime);
#pragma omp declare simd
a5err mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, real* rand5);
#pragma omp declare simd
a5err mccc_wiener_clean(mccc_wienarr* w, real t);

void mccc_wiener_boxmuller(real* randVar, int Ndim);

#pragma omp end declare target

#endif
