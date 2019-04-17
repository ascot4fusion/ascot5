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
enum {
    endcond_tmax   = 0x1,   /**< Maximum simulation time */
    endcond_emin   = 0x2,   /**< Minimum energy          */
    endcond_therm  = 0x4,   /**< Thermalized             */
    endcond_wall   = 0x8,   /**< Wall collision          */
    endcond_rhomin = 0x10,  /**< Minimum rho             */
    endcond_rhomax = 0x20,  /**< Maximum rho             */
    endcond_polmax = 0x40,  /**< Poloidal limit          */
    endcond_tormax = 0x80,  /**< Toroidal limit          */
    endcond_cpumax = 0x100, /**< Wall time exceeded      */
    endcond_hybrid = 0x200  /**< Hybrid mode condition   */
};

#pragma omp declare target
void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i,
                      sim_data* sim);
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i,
                      sim_data* sim);
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i,
                      sim_data* sim);
#pragma omp end declare target

void endcond_parse(int endcond, int* endconds);
void endcond_parse2str(int endcond, char* str);

#endif
