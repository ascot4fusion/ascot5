/**
 * @file diag_transcoef.h
 * @brief Header file for diag_transcoef.c.
 *
 * Contains definitions for transport coefficient data structures.
 */
#ifndef DIAG_TRANSCOEF_H
#define DIAG_TRANSCOEF_H

#include "defines.h"
#include "particle.h"
#include "options.h"
#include "diag.h"


void diag_transcoef_init(diag_transcoef_data* data, sim_parameters* params,
                         size_t nmarkers);
void diag_transcoef_free(diag_transcoef_data* data);
void diag_transcoef_update_fo(diag_transcoef_data* data, sim_parameters* params,
                              particle_simd_fo* p_f, particle_simd_fo* p_i);
void diag_transcoef_update_gc(diag_transcoef_data* data, sim_parameters* params,
                              particle_simd_gc* p_f, particle_simd_gc* p_i);
void diag_transcoef_update_ml(diag_transcoef_data* data, sim_parameters* params,
                              particle_simd_ml* p_f, particle_simd_ml* p_i);
#endif
