/**
 * @file coulomb_collisions.h
 * Tools to integrate Coulomb scattering between test particle and background
 * plasma.
 */
#ifndef COULOMB_COLLISIONS_H
#define COULOMB_COLLISIONS_H

#include "data/bfield.h"
#include "data/marker.h"
#include "data/plasma.h"
#include "defines.h"
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
typedef struct
{
    int nextslot[MCCC_NSLOTS]; /**< Integer array where each element shows
                                    where the next wiener process is located.
                                    Indexing starts from 0 and element points
                                    to itself if it is the last element       */
    real time[MCCC_NSLOTS];    /**< Time instances for different Wiener
                                    processes                                 */
    real wiener[MCCC_NDIM * MCCC_NSLOTS]; /**< Ndim x Nslot array of Wiener
                                               process values */
} mccc_wienarr;

DECLARE_TARGET_SIMD
/**
 * @brief Initializes a struct that stores generated Wiener processes
 *
 * @param w Wiener struct to be initialized
 * @param initime time when a Wiener process begins
 */
void mccc_wiener_initialize(mccc_wienarr *w, real initime);

DECLARE_TARGET_SIMD
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
err_t mccc_wiener_generate(mccc_wienarr *w, real t, int *windex, real *rand5);

DECLARE_TARGET_SIMD
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
err_t mccc_wiener_clean(mccc_wienarr *w, real t);

/**
 * @brief Set collision operator data.
 *
 * @param mdata pointer to collision operator data struct
 * @param include_energy can collisions change marker energy, either 0 or 1
 * @param include_pitch  can collisions change marker pitch, either 0 or 1
 * @param include_gcdiff can collisions change GC position, either 0 or 1
 */
void mccc_init(
    mccc_data *mdata, int include_energy, int include_pitch,
    int include_gcdiff);

/**
 * @brief Integrate collisions for one time-step
 *
 * @param p fo struct
 * @param h time-steps for NSIMD markers
 * @param plasma pointer to plasma data
 * @param mdata pointer collision data struct
 * @param rnd array of normally distributed random numbers used to resolve
 *        collisions. Values for marker i are rnd[i*NSIMD + j]
 */
void mccc_go_euler(
    MarkerGyroOrbit *p, real *h, Plasma *plasma, mccc_data *mdata, real *rnd);

/**
 * @brief Integrate collisions for one time-step
 *
 * @param p gc struct
 * @param h time-steps for NSIMD markers
 * @param bfield pointer to magnetic field
 * @param plasma pointer to plasma data
 * @param mdata pointer to collision data struct
 * @param rnd array of normally distributed random numbers used to resolve
 *        collisions. Values for marker i are rnd[i*NSIMD + j]
 */
void mccc_gc_euler(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Plasma *plasma,
    mccc_data *mdata, real *rnd);

/**
 * @brief Integrate collisions for one time-step
 *
 * @param p pointer to gc simd struct
 * @param hin time-steps for NSIMD markers
 * @param hout suggestions for the next timesteps for NSIMD markers
 * @param tol relative error tolerance
 * @param w array holding wiener structs for NSIMD markers
 * @param bfield pointer to magnetic field data
 * @param plasma pointer to plasma data
 * @param mdata pointer to collision data struct
 * @param rnd array of normally distributed random numbers used to resolve
 *        collisions. Values for marker i are rnd[i*NSIMD + j]
 */
void mccc_gc_milstein(
    MarkerGuidingCenter *p, real *hin, real *hout, real tol, mccc_wienarr *w,
    Bfield *bfield, Plasma *plasma, mccc_data *mdata, real *rnd);

#endif
