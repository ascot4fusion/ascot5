/**
 * @file test_plasma_1d.c
 * @brief Test program for 1D plasma evaluation functions
 */
#include <stdio.h>
#include "ascot5.h"
#include "plasma_1d.h"

int main(int argc, char** argv) {
    plasma_1d_offload_data offload_data;
    real* offload_array;
    plasma_1d_init_offload(&offload_data, &offload_array);

    plasma_1d_data data;
    plasma_1d_init(&data, &offload_data, offload_array);

    int species = 1;

    int i;
    for(i = 0; i < data.n_rho; i++) {
        printf("%le %le\n", data.rho[i], data.dens[species*data.n_rho + i]);
    }

    real rho, dens;
    for(rho = 0.0; rho <= 1.1; rho += 0.005) {
        dens = plasma_1d_eval_dens(rho, species, &data);
        printf("%le %le\n", rho, dens);
    }

    return 0;
}
