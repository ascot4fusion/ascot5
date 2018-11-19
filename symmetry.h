/**
 * @file symmetry.h
 * @brief Header file for symmetry.c
*/
#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "ascot5.h"
#include "error.h"

typedef enum symmetry_type {
    symmetry_type_none, symmetry_type_periodic, symmetry_type_stellarator
} symmetry_type;

#pragma omp declare target
a5err symmetry_apply_scalar(real r_phi_z[], real r0, real phi0, real z0,
                            symmetry_type type, real period_length);
a5err symmetry_apply_vector(real r_phi_z[], real scaling[], real r0, real phi0, real z0,
                            symmetry_type type, real period_length);
#pragma omp end declare target

#endif
