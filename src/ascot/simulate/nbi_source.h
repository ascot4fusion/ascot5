/**
 * @file nbi_source.h
 * Model for the neutral beam injection that calculates shinethrough and birth
 * profile of NBI ions.
 */
#ifndef NBI_SOURCE_H
#define NBI_SOURCE_H

#include "data/marker.h"
#include "data/nbi.h"
#include "datatypes.h"
#include "defines.h"
#include <stddef.h>

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
 * @param inj Injector data.
 * @param tlim Time interval over which the markers are injected.
 *        The marker initial time (at the point of injection) is sampled from
 *        this interval.
 * @param nmrk Number of markers to inject.
 * @param mrk Marker struct to fill.
 */
void nbi_source_inject_markers(
    Simulation *sim, Nbi *inj, real tlim[2], size_t nmrk, State *mrk);

/**
 * Trace a neutral marker until it has ionized or hit wall.
 *
 * This function is for the most part identical to simulate_fo with few
 * exceptions relevant for BBNBI.
 *
 * @param pq Pointer to the marker queue containing the initial neutrals.
 * @param sim Pointer to the simu struct with initialized data.
 */
void nbi_source_trace_markers(MarkerQueue *pq, Simulation *sim);

#endif
