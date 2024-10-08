/**
 * @file atomic.h
 * @brief Header file for atomic.c
 */
#ifndef ATOMIC_H
#define ATOMIC_H

#include "../ascot5.h"
#include "../plasma.h"
#include "../neutral.h"
#include "../particle.h"
#include "../random.h"
#include "../asigma.h"

#ifndef GPU 
#pragma omp declare target
#endif
void atomic_fo(particle_simd_fo* p, real* h,
               plasma_data* p_data, neutral_data* n_data,
               random_data* r_data, asigma_data* asigma_data);
#ifndef GPU 
#pragma omp end declare target
#endif
#endif
