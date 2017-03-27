/**
 * @file test_E.c
 * @brief Test program for electric fields
 *
 * Output is written to standard output, redirect with test_E > output.filename
 * Output is in the format:
 *   r, phi, z, rho, drho/dr, drho/dphi, drho/dz, E_r, E_phi, E_z
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include "E_field.h"
#include "B_field.h"
#include "math.h"

int main(int argc, char** argv) {

    if(argc < 10) {
        printf("Usage: test_E nr rmin rmax nz zmin zmax nphi phimin phimax\n");
        exit(1);
    }

    /* Init background */
    B_field_offload_data offload_Bdata;
    real* offload_array;
    B_field_init_offload(&offload_Bdata, &offload_array);

    B_field_data Bdata;
    B_field_init(&Bdata, &offload_Bdata, offload_array);

    E_field_offload_data offload_Edata;
    E_field_init_offload(&offload_Edata, &offload_array);

    E_field_data Edata;
    E_field_init(&Edata, &offload_Edata, offload_array);

    real E[3];
    real rho[4];

    int n_r = atof(argv[1]);
    real r_min = atof(argv[2]);
    real r_max =atof(argv[3]);
    int n_z =atof(argv[4]);
    real z_min =atof(argv[5]);
    real z_max = atof(argv[6]);
    int n_phi = atof(argv[7]);
    real phi_min = atof(argv[8]) / 180 * math_pi;
    real phi_max = atof(argv[9]) / 180 * math_pi;

    real* r = (real*) malloc(n_r * sizeof(real));
    math_linspace(r, r_min, r_max, n_r);
    real* z = (real*) malloc(n_z * sizeof(real));
    math_linspace(z, z_min, z_max, n_z);
    real* phi = (real*) malloc(n_phi * sizeof(real));
    math_linspace(phi, phi_min, phi_max, n_phi);

    int i, j, k;
    for(i = 0; i < n_r; i++) {
        for(j = 0; j < n_z; j++) {
            for(k = 0; k < n_phi; k++) {
                B_field_eval_rho_drho(rho, r[i], phi[k], z[j], &Bdata);
                E_field_eval_E(E, rho, &Edata);
                printf("%le, %le, %le, ", r[i], phi[k], z[j]);
                printf("%le, %le, %le, %le, ", rho[0], rho[1], rho[2], rho[3]);
                printf("%le, %le, %le\n", E[0], E[1], E[2]);
            }
        }
    }
    
    return 0;
}
