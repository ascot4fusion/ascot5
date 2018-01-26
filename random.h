/**
 * @file random.h
 * @brief Header file for random.c
 */
#ifndef RANDOM_H
#define RANDOM_H

#include <stdlib.h>
#define random_init(seed) srand48(seed)
#define random_uniform() drand48()
#define random_uniform_simd(N, r) random_drand48_uniform_simd(N, r)

#endif
