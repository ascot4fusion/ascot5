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

#include "ascot5.h"
#include "simulate.h"

/**
 * Simulate NBI injection.
 *
 * This function initializes neutrals and traces them until they have ionized or
 * hit the wall.
 *
 * @param sim Pointer to the simulation data structure.
 * @param nprt Number of markers to be injected.
 * @param t1 Time instant when the injector is turned on.
 * @param t2 Time instant when the injector is turned off.
 * @param p Pointer to the marker array which is allocated here.
 */
void bbnbi_simulate(
    sim_data *sim, int nprt, real t1, real t2, particle_state **p);

#endif
