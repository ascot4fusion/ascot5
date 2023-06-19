/**
 * @file atomic.h
 * @brief Header file for atomic.c
 *
 * @todo If possible, define more parameters as uniform!
 */
#ifndef ATOMIC_H
#define ATOMIC_H

#include "ascot5.h"
#include "plasma.h"
#include "neutral.h"
#include "particle.h"
#include "random.h"
#include "asigma.h"

#pragma omp declare target
void atomic_fo(particle_simd_fo* p, real* h,
               plasma_data* p_data, neutral_data* n_data,
               random_data* r_data, asigma_data* asgm_data);
#pragma omp declare simd uniform(asgm_data)
a5err atomic_rates(real* rate_eff_ion, real* rate_eff_rec,
                   int z_1, int a_1, real m_1,
                   const int* z_2, const int* a_2, const real* m_2,
                   asigma_data* asgm_data,
                   int q, real E,
                   int N_pls_spec,
                   real* T, real T_0,
                   real* n, real n_0);
#pragma omp declare simd
a5err atomic_react(int* q, real dt,
                   real rate_eff_ion, real rate_eff_rec,
                   int z_1, real rnd);
#pragma omp end declare target

#endif
