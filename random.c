/**
 * @file random.c
 * @brief Random number generator interface
 */

#if defined(RANDOM_MKL)

#include <mkl_vsl.h>
#include "random.h"

void random_mkl_init(random_data* rdata, int seed) {
    vslNewStream(&rdata->r, VSL_BRNG_SFMT19937, seed);
}

double random_mkl_uniform(random_data* rdata) {
    double r;
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rdata->r, 1, &r, 0.0, 1.0);
    return r;
}

void random_mkl_uniform_simd(random_data* rdata, int n, double* r) {
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rdata->r, n, r, 0.0, 1.0);
}

#elif defined(RANDOM_GSL)

#include <gsl/gsl_rng.h>
#include "random.h"

void random_gsl_init(random_data* rdata, int seed) {
    rdata->r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rdata->r, seed);
}

double random_gsl_uniform(random_data* rdata) {
    return gsl_rng_uniform(rdata->r);
}

void random_gsl_uniform_simd(random_data* rdata, int n, double* r) {
    #pragma omp simd
    for(int i = 0; i < n; i++) {
        r[i] = gsl_rng_uniform(rdata->r);
    }
}




#else /* No RNG lib defined, use drand48 */

#define _XOPEN_SOURCE 500

#include <stdlib.h>
#include "random.h"

void random_drand48_uniform_simd(int n, double* r) {
    #pragma omp simd
    for(int i = 0; i < n; i++) {
        r[i] = drand48();
    }
}

#endif
