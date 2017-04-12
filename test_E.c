/**
 * @file test_E.c
 * @brief Test program for electric fields
 * Output is written to standard output, redirect with test_E > output.filename
 */
#include <stdio.h>
#include <stdlib.h>
#include "E_field.h"
#include "B_field.h"
#include "math.h"

int main(int argc, char** argv) {
    if(argc < 10) {
        printf("Usage: test_E nr rmin rmax nphi phimin phimax nz zmin zmax\n");
        exit(1);
    }

    /* Init magnetic background */
    B_field_offload_data offload_Bdata;
    real* offload_array;
    B_field_init_offload(&offload_Bdata, &offload_array);
    B_field_data Bdata;
    B_field_init(&Bdata, &offload_Bdata, offload_array);
    /* Init electric field */    
    E_field_offload_data offload_Edata;
    E_field_init_offload(&offload_Edata, &offload_array);
    E_field_data Edata;
    E_field_init(&Edata, &offload_Edata, offload_array);

    real E[3];
    real rho[4];

    int n_r = atof(argv[1]);
    real r_min = atof(argv[2]);
    real r_max = atof(argv[3]);
    int n_phi = atof(argv[4]);
    real phi_min = atof(argv[5]) / 180 * math_pi;
    real phi_max = atof(argv[6]) / 180 * math_pi;
    int n_z = atof(argv[7]);
    real z_min = atof(argv[8]);
    real z_max = atof(argv[9]);

    real* r = (real*) malloc(n_r * sizeof(real));
    math_linspace(r, r_min, r_max, n_r);
    real* z = (real*) malloc(n_z * sizeof(real));
    math_linspace(z, z_min, z_max, n_z);
    real* phi = (real*) malloc(n_phi * sizeof(real));
    math_linspace(phi, phi_min, phi_max, n_phi);

    /* Write header specifying grid dimensions */
    printf("%d, %le, %le, ", n_r, r_min, r_max);
    printf("%d, %le, %le, ", n_phi, phi_min, phi_max);
    printf("%d, %le, %le\n", n_z, z_min, z_max);

    int i, j, k;
    for(i = 0; i < n_r; i++) {
        for(j = 0; j < n_z; j++) {
            for(k = 0; k < n_phi; k++) {
                B_field_eval_rho_drho(rho, r[i], phi[k], z[j], &Bdata);
                E_field_eval_E(E, rho, &Edata);
                printf("%le %le %le %le ", rho[0], rho[1], rho[2], rho[3]);
                printf("%le %le %le\n", E[0], E[1], E[2]);
            }
        }
    }
    return 0;
}
