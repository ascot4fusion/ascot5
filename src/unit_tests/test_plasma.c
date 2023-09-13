/**
 * @file test_plasma.c
 * @brief Test program for plasma evaluation functions
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../math.h"
#include "../plasma.h"
#include "../hdf5_interface.h"
#include "../offload.h"

/**
 * Main function for the test program.
 */
int main(int argc, char** argv) {

    if(argc < 11) {
        printf("Usage: test_plasma nr rmin rmax nphi phimin phimax nz zmin zmax fname\n");
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
    int*  wall_int_offload_array;
    real* offload_array;
    int n;
    input_particle* p;

    /* read_options(argc, argv, &sim); */
    strcpy(sim.hdf5_in, "ascot.h5");
    strcpy(sim.hdf5_out, "ascot");

    hdf5_interface_read_input(&sim,
                              hdf5_input_plasma,
                              &B_offload_array,
                              &E_offload_array,
                              &plasma_offload_array,
                              &neutral_offload_array,
                              &wall_offload_array,
                              &wall_int_offload_array,
                              NULL,
                              NULL,
                              &p,
                              &n);
    /* Init plasma */
    offload_package offload_data;
    offload_init_offload(&offload_data, &offload_array);
    offload_pack(&offload_data, &offload_array, plasma_offload_array,
                 sim.plasma_offload_data.offload_array_length);

    plasma_data data;
    plasma_init(&data, &sim.plasma_offload_data, plasma_offload_array);

    int species = 1;

    real rho, dens;
    for(rho = 0.0; rho <= 1.1; rho += 0.005) {
        plasma_eval_dens(&dens, rho, 0, 0, 0, 0, species, &data);
        fprintf(f,"%le %le\n", rho, dens);
    }

    fclose(f);
    return 0;
}
