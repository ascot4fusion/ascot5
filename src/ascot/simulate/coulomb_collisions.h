/**
 * @file coulomb_collisions.h
 * @brief Header file for mccc package
 */
#ifndef COULOMB_COLLISIONS_H
#define COULOMB_COLLISIONS_H

#include "data/bfield.h"
#include "defines.h"
#include "data/marker.h"
#include "data/plasma.h"
#include "utils/random.h"

/**
 * @brief Defines minimum energy boundary condition
 *
 * This times local electron temperature is minimum energy boundary. If guiding
 * center energy goes below this, it is mirrored to prevent collision
 * coefficients from diverging.
 */
#define MCCC_CUTOFF 0.1



/**
 * Wiener process dimension. NDIM=5 because only guiding centers are simulated
 * with adaptive time step
 */
#define MCCC_NDIM 5

/**
 * Maximum slots in Wiener array which means this is the maximum number of time
 * step reductions.
 */
#define MCCC_NSLOTS WIENERSLOTS

/**
 * @brief Struct for storing Wiener processes.
 *
 * Elements of this struct should not be changed outside mccc package.
 */
typedef struct {
    int nextslot[MCCC_NSLOTS]; /**< Integer array where each element shows
                                    where the next wiener process is located.
                                    Indexing starts from 0 and element points
                                    to itself if it is the last element       */
    real time[MCCC_NSLOTS];    /**< Time instances for different Wiener
                                    processes                                 */
    real wiener[MCCC_NDIM*MCCC_NSLOTS]; /**< Ndim x Nslot array of Wiener
                                             process values                   */
} mccc_wienarr;

DECLARE_TARGET_SIMD
void mccc_wiener_initialize(mccc_wienarr* w, real initime);
DECLARE_TARGET_SIMD
err_t mccc_wiener_generate(mccc_wienarr* w, real t, int* windex, real* rand5);
DECLARE_TARGET_SIMD
err_t mccc_wiener_clean(mccc_wienarr* w, real t);


void mccc_init(
    mccc_data *mdata, int include_energy, int include_pitch,
    int include_gcdiff);
void mccc_go_euler(
    MarkerGyroOrbit *p, real *h, Plasma *plasma, mccc_data *mdata,
    real *rnd);
void mccc_gc_euler(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Plasma *plasma,
    mccc_data *mdata, real *rnd);
void mccc_gc_milstein(
    MarkerGuidingCenter *p, real *hin, real *hout, real tol, mccc_wienarr *w,
    Bfield *bfield, Plasma *plasma, mccc_data *mdata, real *rnd);

#endif
