/**
 * @file bbnbi5.h
 * BBNBI5 models neutral beam injectors and is used to evaluate shine-through
 * and beam birth-profile.
 *
 * Neutral markers are generated from injector geometry
 * and traced until they ionize or hit the wall. Several injectors can be
 * modelled simultaneously keeping in mind that in this case the output
 * the injector from which a particle originated is lost.
 */
#ifndef BBNBI5_H
#define BBNBI5_H

#include "defines.h"
#include "simulate.h"

/**
 * Inject neutrals from an injector.
 *
 * This function initializes neutral markers at the beamlet positions and
 * launches them in a (random) direction based on the injector specs.
 * The marker is traced until it enters the magnetic field, at which point
 * the particle struct is filled with only the particle data, and the struct
 * is returned.
 *
 * @param p Pointer where generated markers are stored.
 * @param nprt Number of markers to be injected or generated.
 * @param ngenerated Number of markers that have already been generated.
 * @param t0 Time when the injector is turned on.
 * @param t1 Time when the injector is turned off.
 * @param inj Pointer to injector data.
 * @param sim Pointer to the sim struct with initialized data.
 */
void bbnbi_inject_markers(
    particle_state *p, int nprt, int ngenerated, real t0, real t1,
    nbi_injector *inj, sim_data *sim);

/**
 * Trace a neutral marker until it has ionized or hit wall.
 *
 * This function is for the most part identical to simulate_fo with few
 * exceptions relevant for BBNBI.
 *
 * @param pq Pointer to the marker queue containing the initial neutrals.
 * @param sim Pointer to the simu struct with initialized data.
 */
void bbnbi_trace_markers(particle_queue *pq, sim_data *sim);

#endif
