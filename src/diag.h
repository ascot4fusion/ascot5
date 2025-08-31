/**
 * @file diag.h
 * @brief Header file for diag.c
 */
#ifndef DIAG_H
#define DIAG_H
#include "ascot5.h"
#include "particle.h"
#include "B_field.h"
#include "simulate.h"

int diag_init(sim_data* sim, int Nmrk);
void diag_free(sim_data* sim);
void diag_offload(sim_data* sim);

void diag_update_fo(sim_data* sim, B_field_data* Bdata, particle_simd_fo* p_f,
                    particle_simd_fo* p_i);

void diag_update_gc(sim_data* sim, B_field_data* Bdata, particle_simd_gc* p_f,
                    particle_simd_gc* p_i);

void diag_update_ml(sim_data* sim, particle_simd_ml* p_f,
                    particle_simd_ml* p_i);


#endif
