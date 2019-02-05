/**
 * @file test_E.c
 * @brief Test program for electric fields
 * Output is written to standard output, redirect with test_E > output.filename
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../math.h"
#include "../consts.h"
#include "../B_field.h"
#include "../E_field.h"
#include "../hdf5_interface.h"
#include "../offload.h"

int main(int argc, char** argv) {

    if(argc < 11) {
        printf("Usage: test_E nr rmin rmax nphi phimin phimax nz zmin zmax fname\n");
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

    hdf5_interface_read_input(&sim, &B_offload_array, &E_offload_array,
                              &plasma_offload_array,
                              &neutral_offload_array,
                              &wall_offload_array, &p, &n);

    /* Init magnetic background */
    offload_package offload_data;
    offload_init_offload(&offload_data, &offload_array);
    offload_pack(&offload_data, &offload_array, B_offload_array,
                 sim.B_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, E_offload_array,
                 sim.E_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, plasma_offload_array,
                 sim.plasma_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, wall_offload_array,
                 sim.wall_offload_data.offload_array_length);

    B_field_data Bdata;
    B_field_init(&Bdata, &sim.B_offload_data, B_offload_array);
    E_field_data Edata;
    E_field_init(&Edata, &sim.E_offload_data, E_offload_array);

    real E[3];
    real rho[4];

    int n_r = atof(argv[1]);
    real r_min = atof(argv[2]);
    real r_max = atof(argv[3]);
    int n_phi = atof(argv[4]);
    real phi_min = atof(argv[5]) / 180 * CONST_PI;
    real phi_max = atof(argv[6]) / 180 * CONST_PI;
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
    fprintf(f,"%d %le %le ", n_r, r_min, r_max);
    fprintf(f,"%d %le %le ", n_phi, phi_min, phi_max);
    fprintf(f,"%d %le %le\n", n_z, z_min, z_max);

    int i, j, k;
    for(i = 0; i < n_r; i++) {
        for(j = 0; j < n_phi; j++) {
            for(k = 0; k < n_z; k++) {
                B_field_eval_rho_drho(rho, r[i], phi[j], z[k], &Bdata);
                /* Correct Jacobian */
                rho[2] = rho[2]/r[i];
                E_field_eval_E(E, r[i], phi[j], z[k], 0, &Edata, &Bdata);
                fprintf(f,"%le %le %le %le ", rho[0], rho[1], rho[2], rho[3]);
                fprintf(f,"%le %le %le\n", E[0], E[1], E[2]);
            }
        }
    }
    return 0;
}
