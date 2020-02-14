/**
 * @file test_biosaw.c
 * @brief Test program for biosaw
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../biosaw.h"
#include "../consts.h"

int main(int argc, char** argv) {
    int n = 100;

    real* coil_x = (real*) malloc(n*sizeof(real));
    real* coil_y = (real*) malloc(n*sizeof(real));
    real* coil_z = (real*) malloc(n*sizeof(real));

    for(int i = 0; i < n-1; i++) {
        coil_x[i] = cos(2*i*CONST_PI/(n-1));
        coil_y[i] = 0;
        coil_z[i] = sin(2*i*CONST_PI/(n-1));
    }

    coil_x[n-1] = coil_x[0];
    coil_y[n-1] = coil_y[0];
    coil_z[n-1] = coil_z[0];

    real x[3] = {0};
    real B[3];
    biosaw_calc_B(x, n, coil_x, coil_y, coil_z, B);

    printf("%le %le %le\n", B[0], B[1], B[2]);
}
