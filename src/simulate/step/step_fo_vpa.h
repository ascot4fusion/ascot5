/**
 * @file step_fo_vpa.h
 * @brief Header file for step_fo_vpa.c
 */
#ifndef STEP_FO_VPA_H
#define STEP_FO_VPA_H
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../boozer.h"
#include "../../mhd.h"
#include "../../particle.h"
#include "../../rf_fields_fo.h"

// Defining the a type that describes a generic function to push the particles.
typedef void (*push_fo_fnt)(particle_simd_fo*, real*, B_field_data*, E_field_data*, RF2D_fields*);

// 2nd order integrators, without phase correction.
void step_fo_vpa(particle_simd_fo* p, real* h, B_field_data* Bdata,
                 E_field_data* Edata, RF2D_fields* rffield_data);
void step_fo_vpa_mhd(particle_simd_fo* p, real* h, B_field_data* Bdata,
                     E_field_data* Edata,RF2D_fields* rffield_data, boozer_data* boozer, 
                     mhd_data* mhd);

// With phase correction.
void step_fo_vpa_full(particle_simd_fo* p, real* h, B_field_data* Bdata,
                      E_field_data* Edata, RF2D_fields* rffield_data);
void step_fo_vpa_borisA(particle_simd_fo* p, real* h, B_field_data* Bdata,
                      E_field_data* Edata, RF2D_fields* rffield_data);

// Higher order integrators.
void step_fo_vpa_4th(particle_simd_fo* p, real* h, B_field_data* Bdata,
                     E_field_data* Edata, RF2D_fields* rffield_data);

#endif
