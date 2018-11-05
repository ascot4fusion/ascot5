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
 * @brief Marker end condition masks
 *
 * These masks are used to mark specific end condition as being active.
 */
enum {
    endcond_tmax   = 0b1, /**< Maximum simulation time */
    endcond_emin   = 0b10, /**< Minimum energy */
    endcond_therm  = 0b100, /**< Thermalized */
    endcond_wall   = 0b1000, /**< Wall collision */
    endcond_rhomin = 0b10000, /**< Minimum rho */
    endcond_rhomax = 0b100000, /**< Maximum rho */
    endcond_polmax = 0b1000000, /**< Poloidal limit */
    endcond_tormax = 0b10000000, /**< Toroidal limit */
    endcond_cpumax = 0b100000000, /**< Wall time exceeded */
    endcond_hybrid = 0b1000000000 /**< Hybrid mode condition */
};

/**
 * @brief Human-readable end conditions
 *
 * @deprecated Intended for human-readable output. However, these don't work as
 * multiple end conditions can be active simultaneously.
 *
 * @todo See if this is used anywhere and remove this
 */
enum {
    endcond_id_rejected = -2,
    endcond_id_aborted = -1,
    endcond_id_tmax = 1,
    endcond_id_emin = 2,
    endcond_id_wall = 3,
    endcond_id_therm = 4,
    endcond_id_cpumax = 5,
    endcond_id_rhomin = 10,
    endcond_id_rhomax = 11,
    endcond_id_polmax = 12,
    endcond_id_tormax = 13
};

#pragma omp declare target
void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i,
                      sim_data* sim);
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i,
                      sim_data* sim);
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i,
                      sim_data* sim);
#pragma omp end declare target

#endif
