/**
 * @file nbi_source.h
 * Model for the neutral beam injection that calculates shinethrough and birth
 * profile of NBI ions.
 */
#ifndef NBI_SOURCE_H
#define NBI_SOURCE_H

#include "data/marker.h"
#include "datatypes.h"
#include "defines.h"
#include "simulate.h"
#include <stddef.h>

/**
 * Sample initial phase-space coordinates for an injected marker.
 *
 * @param xyz Sampled positional coordinates in Cartesian basis [m].
 * @param vxyz Sampled velocity coordinates in Cartesian basis [m/s].
 * @param rnd Array of random numbers where first three are uniformly
 *        distributed between (0, 1] and the last two are normally distributed.
 * @param nbi The injector data.
 */
void nbi_source_sample_beam(real xyz[3], real vxyz[3], real rnd[5], Nbi *nbi);

/**
 * Inject neutrals from an injector.
 *
 * This function initializes neutral markers at the beamlet positions. The
 * initial velocity is randomly sampled from a distribution specified by the
 * injector data. Setting marker IDs are left to the caller.
 *
 * The marker is traced until it enters the magnetic field, at which point
 * the marker struct is filled with only the particle data (i.e. guiding center
 * position remains unfilled).
 *
 * @param sim Simulation data.
 * @param nbi Injector data.
 * @param time_step Simulation time-step [s].
 * @param max_dist Maximum distance the marker can travel without entering the
 *        magnetic field before it is aborted.
 * @param tlim Time interval from which markers' initial time is sampled
 *        (separately for each marker).
 * @param nmrk Number of markers to inject.
 * @param mrk Marker structs to fill.
 */
void nbi_source_inject_markers(
    Simulation *sim, Nbi *nbi, real time_step, real max_dist, real tlim[2],
    size_t nmrk, State mrk[nmrk]);

#endif
