/**
 * @file test_random.c
 * @brief Test program for random number generator
 */
#include <stdio.h>
#include <omp.h>
#include "ascot5.h"
#include "random.h"

#define N 100000

int main(int argc, char** argv) {
    double r[N];

    random_init(12345);

    double t1, t2, t3, t4;
    t1 = omp_get_wtime();

    for(int i = 0; i < N; i++) {
        r[i] = random_uniform();
    }

    t2 = omp_get_wtime();

    random_uniform_simd(N, r);

    t3 = omp_get_wtime();

    printf("Serial %lf, SIMD %lf\n", t2-t1, t3-t2);

    return 0;
}
