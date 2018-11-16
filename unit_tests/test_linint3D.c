/**
 * @file test_linint3D.c
 * @brief Test program for trilinear interpolation
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../math.h"
#include "../linint/linint3D.h"

int main(int argc, char** argv) {

    if(argc < 2) {
        printf("Usage: test_linint3D fname\n");
        exit(1);
    }

    FILE* f = fopen(argv[1], "w");
    
    real val;

    int n_r = 14;
    real r_min = 0.0;
    real r_max = 2.0;
    real r_grid = (r_max - r_min)/(n_r-1);
    int n_phi = 4;
    real phi_min = 0.0;
    real phi_max = 2*CONST_PI;
    real phi_grid = (phi_max - phi_min)/(n_phi-1);
    int n_z = 70;
    real z_min = -2.0;
    real z_max = 2.0;
    real z_grid = (z_max - z_min)/(n_z-1);
    
    real* r = (real*) malloc(n_r*10 * sizeof(real));
    math_linspace(r, r_min, r_max, n_r*10);
    real* phi = (real*) malloc(n_phi*10 * sizeof(real));
    math_linspace(phi, phi_min, phi_max, n_phi*10);
    real* z = (real*) malloc(n_z*10 * sizeof(real));
    math_linspace(z, z_min, z_max, n_z*10);
    
    /* Initialize test data */
    real* fn = (real*) malloc(n_r * n_z * n_phi * sizeof(real));
    for (int k = 0; k < n_phi; k++) {
        for (int j = 0; j < n_z; j++) {
            for (int i = 0; i < n_r; i++) {
                /* Constant for phi, varying with r and z */
                fn[k*n_z*n_r + j*n_r + i] = z[j*10]*r[i*10];
            }
        }
    }

    linint3D_data str;
    linint3D_init(&str, fn, n_r, n_phi, n_z,
                  r_min, r_max, r_grid,
                  phi_min, phi_max, phi_grid,
                  z_min, z_max, z_grid);

    for (int k = 0; k < n_phi*10; k++) {
        for (int j = 0; j < n_z*10; j++) {
            for (int i = 0; i < n_r*10; i++) {
                linint3D_eval(&val, &str, r[i], phi[k], z[j]); 
                fprintf(f,"%le %le %le %le\n", r[i], z[j], phi[k], val);
            }
        }
    }

    fclose(f);
    return 0;
}
