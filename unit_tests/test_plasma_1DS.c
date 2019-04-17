/**
 * @file test_plasma_1DS.c
 * @brief Test program for 1D plasma evaluation functions
 */
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../plasma_1DS.h"
#include "../hdf5io/hdf5_plasma.h"
#include "../hdf5.h"
#include "../hdf5io/hdf5_helpers.h"
#include "../hdf5_hl.h"
#include "../math.h"

int main(int argc, char** argv) {
    plasma_1DS_offload_data offload_data;
    real* offload_array;
    hid_t f = hdf5_open("ascot.h5");
    hdf5_plasma_init_offload_1DS(f, &offload_data, &offload_array);

    plasma_1DS_data data;
    plasma_1DS_init(&data, &offload_data, offload_array);

    int nrho = 1000;
    real* rho = (real*) malloc(nrho * sizeof(real));
    math_linspace(rho, 0.0, 1.0, nrho);
    real denss[data.n_species];
    real temps[data.n_species];
    for (int i = 0; i < nrho; i++) {
        printf("%le ", rho[i]);
        plasma_1DS_eval_densandtemp(denss, temps, rho[i], &data);
        for (int species = 0; species < data.n_species; species++) {
            printf("%le %le ", temps[species], denss[species]);
        }
        printf("\n");
    }
    plasma_1DS_free(&data);

    return 0;
}
