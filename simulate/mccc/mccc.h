/**
 * @file mccc.h
 * @brief Header file for mccc.c
*/
#ifndef MCCC_H
#define MCCC_H

#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../particle.h"
#include "../../random.h"
#include "mccc_wiener.h"

#pragma omp declare target
void mccc_setoperator(int include_energy, int include_pitch,
                      int include_gcdiff);
void mccc_update_fo(particle_simd_fo* p, B_field_data* Bdata, plasma_data* pdata, real* coldata, 
                    real* clogab, real* F, real* Dpara, real* Dperp, real* K, real* nu);
void mccc_update_gc(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, real* coldata,
                    real* clogab, real* Dpara, real* DX, real* K, real* nu, real* dQ, real* dDpara);
void mccc_collfreq_gc(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, real* coldata, 
                      real* nu, int i);
void mccc_step_fo_fixed(particle_simd_fo* p, B_field_data* Bdata, plasma_data* pdata, random_data* rdata, real* coldata, real* h);
void mccc_step_gc_fixed(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, random_data* rdata, real* coldata, real* h);
void mccc_step_gc_adaptive(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, random_data* rdata, real* coldata, real* hin, real* hout, mccc_wienarr** w, real tol);
void mccc_printerror(int err);
#pragma omp end declare target

#endif
