/**
 * @file endcond.h
 * @brief Header file for endcond.c
 */
#ifndef ENDCOND_H
#define ENDCOND_H

#include "particle.h"
#include "simulate.h"

enum {
    endcond_tmax   = 0b1,
    endcond_emin   = 0b10,
    endcond_therm  = 0b100,
    endcond_wall   = 0b1000,
    endcond_rhomax = 0b10000
};

enum {
    endcond_id_tmax = 1,
    endcond_id_emin = 2,
    endcond_id_wall = 3,
    endcond_id_therm = 4,
    endcond_id_maxrho = 10
};

#pragma omp declare target
void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i, sim_data* sim);
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i, sim_data* sim);
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i, sim_data* sim);
#pragma omp end declare target

#endif
