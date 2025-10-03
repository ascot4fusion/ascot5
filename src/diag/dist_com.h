/**
 * @file dist_com.h
 * @brief Header file for dist_com.c
 */
#ifndef DIST_COM_H
#define DIST_COM_H

#include <stdlib.h>
#include "ascot5.h"
#include "particle.h"
#include "B_field.h"
#include "diag.h"



int dist_COM_init(dist_COM_data* data);
void dist_COM_free(dist_COM_data* data);
void dist_COM_offload(dist_COM_data* data);
void dist_COM_update_fo(dist_COM_data* dist, B_field_data*Bdata,
                        particle_simd_fo* p_f, particle_simd_fo* p_i);
void dist_COM_update_gc(dist_COM_data* dist, B_field_data* Bdata,
                        particle_simd_gc* p_f, particle_simd_gc* p_i);

#endif
