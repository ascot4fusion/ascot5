/**
 * @file step_gc_rk4.h
 * @brief Header file for step_gc_rk4.c
 */
#ifndef STEP_GC_CASHKARP_H
#define STEP_GC_CASHKARP_H

#include "../../B_field.h"
#include "../../E_field.h"
#include "../../boozer.h"
#include "../../mhd.h"
#include "../../particle.h"

void step_gc_cashkarp(particle_simd_gc* p, real* h, real* hnext, real tol,
                      B_field_data* Bdata, E_field_data* Edata);
void step_gc_cashkarp_mhd(particle_simd_gc* p, real* h, real* hnext, real tol,
                          B_field_data* Bdata, E_field_data* Edata,
                          boozer_data* boozer, mhd_data* mhd);

#endif
