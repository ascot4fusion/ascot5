/**
 * @file step_ml_cashkarp.h
 * @brief Header file for step_ml_cashkarp.c
 */
#ifndef STEP_ML_CASHKARP_H
#define STEP_ML_CASHKARP_H

#include "B_field.h"
#include "particle.h"

#pragma omp declare target
void step_ml_cashkarp(particle_simd_ml* p, real* h, real* hnext, 
		      real tol, B_field_data* Bdata);
#pragma omp end declare target

#endif
