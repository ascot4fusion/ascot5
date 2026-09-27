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
 * Defines minimum energy boundary condition
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
 * Struct for storing Wiener processes.
 *
 * Elements of this struct should not be changed outside mccc package.
 */
typedef struct
{
    /**
     * Integer array where each element shows where the next wiener process is
     * located.
     *
     * Indexing starts from 0 and element pointsto itself if it is the last
     * element.
     */
    int nextslot[MCCC_NSLOTS];

    /**
     * Time instances for different Wiener processes.
     */
    real time[MCCC_NSLOTS];

    /**
     * Ndim x Nslot array of Wiener process values.
     */
    real wiener[MCCC_NDIM * MCCC_NSLOTS];
} mccc_wienarr;

DECLARE_TARGET_SIMD
/**
 * Initialize a struct that stores generated Wiener processes.
 *
 * @param w Wiener struct to be initialized.
 * @param initime Time when a Wiener process begins.
 */
void mccc_wiener_initialize(mccc_wienarr *w, real initime);

/**
 * Offload a struct that stores generated Wiener processes.
 *
 * @param w Wiener struct to be offloaded.
 * @param vector_size The simulation vector size.
 */
void mccc_wiener_offload(mccc_wienarr *w, size_t vector_size);

GPU_DECLARE_TARGET_SIMD
/**
 * Generates a new Wiener process at a given time instant
 *
 * Generates a new Wiener process. The generated process is drawn from
 * normal distribution unless there exists a Wiener process at future
 * time-instance, in which case the process is created using the Brownian
 * bridge.
 *
 * @param w Array that stores the Wiener processes.
 * @param t Time for which the new process will be generated.
 * @param windex Index of the generated Wiener process in the Wiener array.
 * @param rand5 Array of 5 normal distributed random numbers.
 *
 * @return Zero if generation succeeded.
 */
err_t mccc_wiener_generate(mccc_wienarr *w, real t, int *windex, real *rand5);

GPU_DECLARE_TARGET_SIMD
/**
 * Removes Wiener processes from the array that are no longer required.
 *
 * Processes W(t') are redundant if t' <  t, where t is the current simulation
 * time. Note that W(t) should exist before W(t') are removed. This routine
 * should be called each time when simulation time is advanced.
 *
 * @param w Array that stores the Wiener processes.
 * @param t Time for which the new process will be generated.
 *
 * @return Zero if cleaning succeeded.
 */
err_t mccc_wiener_clean(mccc_wienarr *w, real t);

/**
 * Set collision data.
 *
 * @param mdata Collision data.
 * @param include_energy Toggle whether collisions change marker energy.
 * @param include_pitch Toggle whether collisions change marker pitch.
 * @param include_gcdiff Toggle whether collisions change guiding-center
 *        position.
 */
void mccc_init(
    mccc_data *mdata, int include_energy, int include_pitch,
    int include_gcdiff);

/**
 * Integrate collisions for one time-step
 *
 * @param p Marker struct.
 * @param h Time step.
 * @param plasma Plasma data.
 * @param mdata Collision data.
 * @param rnd Array of normally distributed random numbers used to resolve
 *        collisions.
 *
 *        Values for marker i are rnd[i*NSIMD + j]
 */
void mccc_go_euler(
    MarkerGyroOrbit *p, real *h, Plasma *plasma, mccc_data *mdata, real *rnd);

/**
 * Integrate collisions for one time-step
 *
 * @param p Marker struct.
 * @param h Time step.
 * @param bfield Magnetic field data.
 * @param plasma Plasma data.
 * @param mdata Collision data.
 * @param rnd Array of normally distributed random numbers used to resolve
 *        collisions.
 *
 *        Values for marker i are rnd[i*NSIMD + j].
 */
void mccc_gc_euler(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Plasma *plasma,
    mccc_data *mdata, real *rnd);

/**
 * Integrate collisions for one time-step
 *
 * @param p Marker struct.
 * @param hin Time-step.
 * @param acc Acceleration data.
 * @param collfreq Collision frequency to be stored.
 * @param hout Suggestion for the next timestep.
 * @param tol Relative error tolerance
 * @param w Array holding wiener processes.
 * @param bfield Magnetic field data.
 * @param plasma Plasma data.
 * @param mdata Collision data.
 * @param rnd Array of normally distributed random numbers used to resolve
 *        collisions.
 *
 *        Values for marker i are rnd[i*NSIMD + j].
 */
void mccc_gc_milstein(
    MarkerGuidingCenter *p, real *hin, real *acc, real *collfreq, real *hout,
    real tol, mccc_wienarr *w, Bfield *bfield, Plasma *plasma, mccc_data *mdata,
    real *rnd);

#endif
