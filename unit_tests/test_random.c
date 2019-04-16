/**
 * @file test_random.c
 * @brief Test program for random number generator
 */
#include <stdio.h>
#include <omp.h>
#include "../ascot5.h"
#include "../random.h"

#define N 1000000 /**< Number of random numbers to be genrated */

/**
 * Main function for the test program
 */
int main(int argc, char** argv) {
    random_data rdata;

    double r[N];

    random_init(&rdata, 12345);

    double t1, t2, t3;
    t1 = omp_get_wtime();

    for(int i = 0; i < N; i++) {
        r[i] = random_uniform(&rdata);
    }

    t2 = omp_get_wtime();

    random_uniform_simd(&rdata, N, r);

    t3 = omp_get_wtime();

    printf("Serial %lf, SIMD %lf\n", t2-t1, t3-t2);

/*    for(int i = 0; i < N; i++) {
        printf("%le\n", r[i]);
    }*/

    return 0;
}
