/**
 * @file test_hdf5.c
 * @brief Test program for HDF5 functions
 *
 * This test program writes the intial distribution or state of test particles
 * in HDF5 format.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "distributions.h"
#include "B_field.h"
#include "ascot4_interface.h"
#include "hdf5_particlestate.h"
#include "hdf5_helpers.h"
#include "particle.h"

int main(int argc, char** argv) {
    int n;
    particle* p;
    ascot4_read_particles(&p, &n, "input.particles");

    hid_t f = hdf5_create("test.h5");
    hdf5_particlestate_write(f, "inistate", n, p);

/*
    B_field_offload_data B_offload_data;
    real* B_offload_array;
    B_field_init_offload(&B_offload_data, &B_offload_array);

    B_field_data B_data;
    B_field_init(&B_data, &B_offload_data, B_offload_array);

    dist_rzvv_offload_data dist_offload_data;
    real* dist_offload_array;
    dist_rzvv_init_offload(&dist_offload_data, &dist_offload_array,
                           25, 3, 8.5, 35, -4.25, 3.6,
                           25, -1.5e7, 1.5e7, 35, 0, 1.5e7);

    dist_rzvv_data dist_data;
    dist_rzvv_init(&dist_data, &dist_offload_data, dist_offload_array);

    int i, j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < NSIMD; j++) {
            real B[3][NSIMD];
            B_field_eval_B(j, B, p[i].r[j], p[i].phi[j], p[i].z[j], &B_data);
            p[i].B_r[j] = B[0][j];
            p[i].B_phi[j] = B[1][j];
            p[i].B_z[j] = B[2][j];
        }
        dist_rzvv_update(&dist_data, &p[i], 1);
    }

    ascot4_write_dist_rzvv(&dist_offload_data, dist_offload_array, "test.h5");

//    dist_rzvv_print_rz(&dist_offload_data, dist_offload_array);
//    dist_rzvv_print_vv(&dist_offload_data, dist_offload_array);

    B_field_free_offload(&B_offload_data, &B_offload_array);
    dist_rzvv_free_offload(&dist_offload_data, &dist_offload_array);
*/

    free(p);

    return 0;
}

