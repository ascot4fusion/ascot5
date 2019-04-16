/**
 * @file test_interp1Dcomp.c
 * @brief Test program for 1D spline interpolation
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../math.h"
#include "../spline/interp.h"

/**
 * Main function for the test program.
 */
int main(int argc, char** argv) {

    if(argc < 2) {
        printf("Usage: test_interp1Dcomp fname\n");
        exit(1);
    }

    FILE* f = fopen(argv[1], "w");

    real B = 0;
    real B_dB[3];

    int n_r = 100;
    real r_min = 0.0;
    real r_max = 3.0;
    real r_grid = (r_max - r_min)/(n_r-1);

    real* r = (real*) malloc(n_r * sizeof(real));
    math_linspace(r, r_min, r_max, n_r);

    /* Initialize test data */
    real* fn = (real*) malloc(n_r * sizeof(real));
    for (int i = 0; i < n_r; i++) {
        fn[i] = r[i]*r[i];
    }
    /* Natural boundary conditions */
    int bc_r = 0;
    real* c = malloc(n_r * NSIZE_COMP1D * sizeof(real));
    interp1Dcomp_init_coeff(c, fn, n_r, bc_r, r_min, r_max);
    interp1D_data str;
    interp1Dcomp_init_spline(&str, c, n_r, bc_r, r_min, r_max);

    for (int i = 0; i < n_r; i++) {
        interp1Dcomp_eval_f(&B, &str, r[i]);
        interp1Dcomp_eval_df(B_dB, &str, r[i]);
        fprintf(f,"%le %le ", r[i], r[i]*r[i]);
        /* fprintf(f,"%le %le \n", str.c[i*2], str.c[i*2+1]); */
        fprintf(f,"%le %le %le %le\n", B, B_dB[0], B_dB[1], B_dB[2]);
    }
    fclose(f);
    return 0;
}
