/**
 * @file wall.c
 * @brief Random number generator interface
 */
#include <stdlib.h>
#include "random.h"

void random_drand48_uniform_simd(int n, double* r) {
    #pragma omp simd
    for(int i = 0; i < n; i++) {
        r[i] = drand48();
    }
}
