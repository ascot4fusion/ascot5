/**
 * @file dist_6D.h
 * @brief Header file for dist_6D.c
 */
#ifndef DIST_6D_H
#define DIST_6D_H

#include <stdlib.h>
#include "ascot5.h"
#include "particle.h"
#include "diag.h"



int dist_6D_init(dist_6D_data* data);
void dist_6D_free(dist_6D_data* data);
void dist_6D_offload(dist_6D_data* data);
void dist_6D_update_fo(dist_6D_data* dist,
                       particle_simd_fo* p_f, particle_simd_fo* p_i);
void dist_6D_update_gc(dist_6D_data* dist,
                       particle_simd_gc* p_f, particle_simd_gc* p_i);

#endif
