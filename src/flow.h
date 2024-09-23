/**
 * @file flow.h
 * @brief Header file for flow.c
 */
#ifndef FLOW_H
#define FLOW_H

#include "ascot5.h"
#include "B_field.h"
#include "particle.h"
#include "plasma.h"

void flow_gc_to_plasma(particle_simd_gc* p, B_field_data* Bdata,
                       plasma_data* pdata);
void flow_gc_to_lab(particle_simd_gc* p, B_field_data* Bdata,
                    plasma_data* pdata);
void flow_fo_to_plasma(particle_simd_fo* p, B_field_data* Bdata,
                       plasma_data* pdata);
void flow_fo_to_lab(particle_simd_fo* p, B_field_data* Bdata,
                    plasma_data* pdata);
#endif
