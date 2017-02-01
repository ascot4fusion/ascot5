/**
 * @file interact.h
 * @brief Header file for interact.c
 */
#ifndef INTERACT_H
#define INTERACT_H
#include "B_field.h"
#include "plasma_1d.h"
#include "particle.h"

#pragma omp declare target
void interact_step_fo_euler(particle_simd_fo* p, real t, real h,
                            B_field_data* Bdata, plasma_1d_data* pdata);
void interact_step_gc_euler(particle_simd_gc* p, real t, real h,
                            B_field_data* Bdata, plasma_1d_data* pdata);
void interact_step_gc_euler_ascot4(particle_simd_gc* p, real t, real h,
                            B_field_data* Bdata, plasma_1d_data* pdata);
#pragma omp declare simd
real interact_D_par(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species);
#pragma omp declare simd
void interact_D_par_dD_par(real* D_par, real* dD_par,
                    real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species);
#pragma omp declare simd
real interact_D_perp(real v, real m, real q, real temp[MAX_SPECIES],
                     real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                     real mass[MAX_SPECIES], int n_species);
#pragma omp declare simd
real interact_nu_fo(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species);
#pragma omp declare simd
real interact_nu_gc(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species);
#pragma omp declare simd
real interact_coulomb_log(real v, real m, real q, real temp_s, real dens_s,
                          real charge_s, real mass_s, real debye_length2);
#pragma omp declare simd
real interact_dcoulomb_log(real v, real m, real q, real temp_s, real dens_s,
                          real charge_s, real mass_s, real debye_length2);
#pragma omp declare simd
real interact_debye_length2(real temp[MAX_SPECIES], real dens[MAX_SPECIES],
                            real charge[MAX_SPECIES], int n_species);
#pragma omp declare simd
real interact_psi(real x);
#pragma omp declare simd
real interact_dpsi(real x);
#pragma omp declare simd
real interact_nu_ascot4(real x);
#pragma omp end declare target

#endif
