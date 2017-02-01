/**
 * @file test_B.c
 * @brief Test program for magnetic fields
 */
#include <stdio.h>
#include <stdlib.h>
#include "B_field.h"
#include "math.h"

int main(int argc, char** argv) {
    if(argc < 10) {
        printf("Usage: test_B nr rmin rmax nz zmin zmax nphi phimin phimax\n");
        exit(1);
    }

    B_field_offload_data offload_data;
    real* offload_array;
    B_field_init_offload(&offload_data, &offload_array);

    B_field_data Bdata;
    B_field_init(&Bdata, &offload_data, offload_array);

    real B[3];
    real B_dB[12];
    real psi[1];
    real rho[1];

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
                B_field_eval_B(B, r[i], phi[k], z[j], &Bdata);
                B_field_eval_B_dB(B_dB, r[i], phi[k], z[j], &Bdata);
                B_field_eval_psi(psi, r[i], phi[k], z[j], &Bdata);
                B_field_eval_rho(rho, psi[0], &Bdata);
                printf("%le, %le, %le, ", r[i], phi[k], z[j]);
                printf("%le, %le, %le, ", B[0], B[1], B[2]);
                printf("%le, %le, %le, ", B_dB[0], B_dB[4], B_dB[8]);
                printf("%le, %le, %le, ", B_dB[1], B_dB[2], B_dB[3]);
                printf("%le, %le, %le, ", B_dB[5], B_dB[6], B_dB[7]);
                printf("%le, %le, %le, ", B_dB[9], B_dB[10], B_dB[11]);
                printf("%le, %le\n", psi[0], rho[0]);
            }
        }
    }

    return 0;
}
