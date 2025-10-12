/**
 * @file dist_5D.h
 * @brief Header file for dist_5D.c
 */
#ifndef DIST_5D_H
#define DIST_5D_H

#include <stdlib.h>
#include "defines.h"
#include "marker.h"
#include "options.h"
#include "diag.h"


size_t dist_5D_index(int i_r, int i_phi, int i_z, int i_ppara, int i_pperp,
                     int i_time, int i_q, size_t step_6, size_t step_5,
                     size_t step_4, size_t step_3, size_t step_2,
                     size_t step_1);
int dist_5D_init(dist_5D_data* data, Options* options);
void dist_5D_free(dist_5D_data* data);
void dist_5D_offload(dist_5D_data* data);
void dist_5D_update_fo(dist_5D_data* dist, Options* options,
                       MarkerGyroOrbit* p_f, MarkerGyroOrbit* p_i);
void dist_5D_update_gc(dist_5D_data* dist, Options* options,
                       MarkerGuidingCenter* p_f, MarkerGuidingCenter* p_i);

#endif
