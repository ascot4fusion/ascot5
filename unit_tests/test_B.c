/**
 * @file test_B.c
 * @brief Test program for magnetic fields
 * Output is written to standard output, redirect with test_B > output.filename
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../math.h"
#include "../consts.h"
#include "../B_field.h"
#include "../hdf5_interface.h"
#include "../offload.h"

/**
 * Main function for the magnetic field test program.
 */
int main(int argc, char** argv) {

    if(argc < 11) {
        printf("Usage: test_B nr rmin rmax nphi phimin phimax nz zmin zmax fname\n");
        exit(1);
    }

    FILE* f = fopen(argv[10], "w");

    sim_offload_data sim;
    sim.mpi_rank = 0;
    sim.mpi_size = 1;

    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* neutral_offload_array;
    real* wall_offload_array;
    real* offload_array;
    int n;
    input_particle* p;

    /* read_options(argc, argv, &sim); */
    strcpy(sim.hdf5_in, "ascot.h5");
    strcpy(sim.hdf5_out, "ascot");

    hdf5_interface_read_input(&sim, hdf5_input_bfield,
                              &B_offload_array, &E_offload_array,
                              &plasma_offload_array,
                              &neutral_offload_array,
                              &wall_offload_array, NULL, NULL,  &p, &n);

    /* Init magnetic background */
    offload_package offload_data;
    offload_init_offload(&offload_data, &offload_array);
    offload_pack(&offload_data, &offload_array, B_offload_array,
                 sim.B_offload_data.offload_array_length);

    B_field_data Bdata;
    B_field_init(&Bdata, &sim.B_offload_data, B_offload_array);

    real B[3];
    real B_dB[12];
    real psi;
    real rho;

    int n_r = atof(argv[1]);
    real r_min = atof(argv[2]);
    real r_max = atof(argv[3]);
    int n_phi = atof(argv[4]);
    real phi_min = atof(argv[5]) / 180 * CONST_PI;
    real phi_max = atof(argv[6]) / 180 * CONST_PI;
    int n_z = atof(argv[7]);
    real z_min =atof(argv[8]);
    real z_max = atof(argv[9]);

    real* r = (real*) malloc(n_r * sizeof(real));
    math_linspace(r, r_min, r_max, n_r);
    real* z = (real*) malloc(n_z * sizeof(real));
    math_linspace(z, z_min, z_max, n_z);
    real* phi = (real*) malloc(n_phi * sizeof(real));
    math_linspace(phi, phi_min, phi_max, n_phi);

    /* Write header specifying grid dimensions */
    fprintf(f,"%d %le %le ", n_r, r_min, r_max);
    fprintf(f,"%d %le %le ", n_phi, phi_min, phi_max);
    fprintf(f,"%d %le %le\n", n_z, z_min, z_max);

    int i, j, k;
    real time = 0;
    for(i = 0; i < n_r; i++) {
        for(j = 0; j < n_phi; j++) {
            for(k = 0; k < n_z; k++) {
                B_field_eval_B(B, r[i], phi[j], z[k], time, &Bdata);
                B_field_eval_B_dB(B_dB, r[i], phi[j], z[k], time, &Bdata);
                B_field_eval_psi(&psi, r[i], phi[j], z[k], time, &Bdata);
                B_field_eval_rho(&rho, psi, &Bdata);
                fprintf(f,"%le %le %le ", B[0], B[1], B[2]);
                fprintf(f,"%le %le %le ", B_dB[1], B_dB[2], B_dB[3]);
                fprintf(f,"%le %le %le ", B_dB[5], B_dB[6], B_dB[7]);
                fprintf(f,"%le %le %le ", B_dB[9], B_dB[10], B_dB[11]);
                fprintf(f,"%le\n", rho);
            }
        }
    }
    fclose(f);
    return 0;
}
