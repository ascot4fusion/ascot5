/**
 * @file endcond.h
 * Marker simulation end conditions
 *
 * In the absence of errors, marker simulation is ended when marker meets even
 * one of the active end conditions. User can choose which end conditions are
 * active.
 *
 * The end conditions are:
 * - tlim: Marker time has passed the simulation time limit
 *
 * - emin: Marker energy is below minimum value
 *
 * - therm: Energy is below value derived from local thermal electron energy
 *
 * - wall: Marker has intersected wall
 *
 * - rhomin: Marker has reached minimum rho (normalized poloidal flux) value
 *
 * - rhomax: Marker has reached maximum rho value
 *
 * - polmax: The total cumulative distance marker has travelled poloidally
 *   exceeds maximum value
 *
 * - tormax: The total cumulative distance marker has travelled toroidally
 *   exceeds maximum value
 *
 * - cpumax: Marker simulation has exceeded maximum wall time
 *
 * - neutr: Marker has been neutralized by an atomic reaction
 *
 * - ioniz: Marker has been ionized by an atomic reaction
 *
 * - hybrid: Not an end condition per se but used to notate that the guiding
 *   center simulation will be resumed as a gyro-orbit simulation
 *
 * As magnetic field lines have no energy, emin and therm are never checked for
 * them. Guiding centers are the only markers for which hybrid is checked.
 *
 * In the code, the end conditions are represented as bit arrays with each bit
 * corresponding to a specific end condition. Each marker has a field "endcond",
 * and when marker meets an end condition, the corresponding bit is flagged.
 * This way if marker simultaneously meets several end conditions, all can be
 * flagged.
 *
 * Additionally, when marker meets an end condition, its running state is set to
 * False which notates its simulation should be discontinued. If the end
 * condition is wall collision, the ID of the wall element the marker collided
 * with is stored in the marker fields.
 *
 * @todo Error checking would be a good idea
 */
#ifndef ENDCOND_H
#define ENDCOND_H

#include "data/marker.h"
#include "datatypes.h"

/**
 * Check end conditions for GC markers.
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f Pointer to SIMD struct storing marker states at the end of current
 *        time-step.
 * @param p_i Pointer to SIMD struct storing marker states at the beginning of
 *        current time-step.
 * @param sim Pointer to simulation data struct.
 */
void endcond_check_go(
    MarkerGyroOrbit *p_f, MarkerGyroOrbit *p_i, Simulation *sim);

/**
 * Check end conditions for FO markers.
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f Pointer to SIMD struct storing marker states at the end of current
 *        time-step.
 * @param p_i Pointer to SIMD struct storing marker states at the beginning of
 *        current time-step.
 * @param sim Pointer to simulation data struct.
 */
void endcond_check_gc(
    MarkerGuidingCenter *p_f, MarkerGuidingCenter *p_i, Simulation *sim);

/**
 * Check end conditions for ML markers.
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f Pointer to SIMD struct storing marker states at the end of current
 *        time-step.
 * @param p_i Pointer to SIMD struct storing marker states at the beginning of
 *        current time-step.
 * @param sim Pointer to simulation data struct.
 */
void endcond_check_fl(
    MarkerFieldLine *p_f, MarkerFieldLine *p_i, Simulation *sim);

#endif
