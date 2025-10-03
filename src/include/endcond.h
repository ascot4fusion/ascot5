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

#include "particle.h"
#include "simulate.h"

/**
 * Marker end condition bit masks.
 *
 * These bit masks are used to mark specific end condition as being active.
 */
extern const unsigned int endcond_tlim;
extern const unsigned int endcond_emin;
extern const unsigned int endcond_therm;
extern const unsigned int endcond_wall;
extern const unsigned int endcond_rhomin;
extern const unsigned int endcond_rhomax;
extern const unsigned int endcond_polmax;
extern const unsigned int endcond_tormax;
extern const unsigned int endcond_cpumax;
extern const unsigned int endcond_hybrid;
extern const unsigned int endcond_neutr;
extern const unsigned int endcond_ioniz;

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
    particle_simd_gc *p_f, particle_simd_gc *p_i, sim_data *sim);

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
void endcond_check_fo(
    particle_simd_fo *p_f, particle_simd_fo *p_i, sim_data *sim);

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
void endcond_check_ml(
    particle_simd_ml *p_f, particle_simd_ml *p_i, sim_data *sim);

#endif
