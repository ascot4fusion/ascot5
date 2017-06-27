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
    endcond_rhomin = 0b10000,
    endcond_rhomax = 0b100000,
    endcond_polmax = 0b1000000,
    endcond_tormax = 0b10000000,
    endcond_cpumax = 0b100000000,
    endcond_hybrid = 0b1000000000
};

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
void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i, sim_data* sim);
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i, sim_data* sim);
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i, sim_data* sim);
#pragma omp end declare target

#endif
