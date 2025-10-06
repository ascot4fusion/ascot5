/**
 * @file dist_rho6D.h
 * @brief Header file for dist_rho6D.c
 */
#ifndef DIST_RHO6D_H
#define DIST_RHO6D_H

#include <stdlib.h>
#include "defines.h"
#include "particle.h"
#include "diag.h"



int dist_rho6D_init(dist_rho6D_data* dist_data);
void dist_rho6D_free(dist_rho6D_data* dist_data);
void dist_rho6D_offload(dist_rho6D_data* dist_data);
void dist_rho6D_update_fo(dist_rho6D_data* dist, particle_simd_fo* p_f,
                          particle_simd_fo* p_i);
void dist_rho6D_update_gc(dist_rho6D_data* dist, particle_simd_gc* p_f,
                          particle_simd_gc* p_i);

#endif
