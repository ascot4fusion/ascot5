/**
 * @file test_ascot4_interface.c
 * @brief Test program for ASCOT4 data export functions
 *
 * This program can be used to call ascot4_write_B for exporting a magnetic
 * field to be processed into ASCOT4 input files with prepare_magn_bkg.m.
 */
#include <stdio.h>
#include "ascot5.h"
#include "ascot4_interface.h"
#include "particle.h"

int main(int argc, char** argv) {

    B_field_offload_data offload_data;
    real* offload_array;
    B_field_init_offload(&offload_data, &offload_array);
    B_field_data data;
    B_field_init(&data, &offload_data, offload_array);

    ascot4_write_B(&data);

/*
    particle_simd* p;
    int n;
    ascot4_read_particles(&p, &n, "input.particles.16");
    int i, j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < NSIMD; j++) {
            printf("%ld %le %le %le %le %le %le %le %le\n",
                   p[i].id[j],
                   p[i].r[j], p[i].phi[j],
                   p[i].z[j], p[i].rdot[j],
                   p[i].phidot[j], p[i].zdot[j],
                   p[i].mass[j], p[i].charge[j]);
        }
    }
*/

    return 0;
}
