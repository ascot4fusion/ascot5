/**
 * @file mccc.h
 * @brief Header file for mccc package
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
void mccc_eval_coefs(real m, real q, real r, real phi, real z, real t, real v,
                     int nv, plasma_data* pdata, real* F, real* Dpara,
                     real* Dperp, real* K, real* nu, int* err);
void mccc_fo_euler(particle_simd_fo* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, real* coldata);
void mccc_gc_euler(particle_simd_gc* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, real* coldata);
void mccc_gc_milstein(particle_simd_gc* p, real* hin, real* hout, real tol,
                      mccc_wienarr** wienarr, B_field_data* Bdata,
                      plasma_data* pdata, random_data* rdata, real* coldata);
#pragma omp end declare target

#endif
