/**
 * @file endcond.h
 * @brief Header file for endcond.c
 *
 * Contains a list declaring all end conditions.
 */
#ifndef ENDCOND_H
#define ENDCOND_H

#include "particle.h"
#include "simulate.h"

/**
 * @brief Marker end condition bit masks
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


void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i,
                      sim_data* sim);
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i,
                      sim_data* sim);
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i,
                      sim_data* sim);
#endif
