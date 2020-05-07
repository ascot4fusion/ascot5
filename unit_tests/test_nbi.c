/**
 * @file test_nbi.c
 * @brief Test program for NBI functions
 */
#define _XOPEN_SOURCE 500 /**< drand48 requires POSIX 1995 standard */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ascot5.h"
#include "../B_field.h"
#include "../hdf5_interface.h"
#include "../hdf5io/hdf5_nbi.h"
#include "../hdf5io/hdf5_helpers.h"
#include "../hdf5io/hdf5_bfield.h"
#include "../hdf5io/hdf5_plasma.h"
#include "../random.h"
#include "../nbi.h"
#include "../particle.h"
#include "../plasma.h"
#include "../suzuki.h"

int main(int argc, char** argv) {
/*    int Z[] = {1,4};
    int A[] = {3,9};
    real ni[] = {5e19, 1e17};
    printf("%le\n", suzuki_sigmav(150, 5.0025e19, 5000, 2, ni, A, Z));

    return 0;*/

    sim_offload_data sim;
    strcpy(sim.hdf5_in, "ascot.h5");

    real* B_offload_array;
    real* plasma_offload_array;
    real* wall_offload_array;
    hdf5_interface_read_input(&sim, hdf5_input_bfield | hdf5_input_plasma |
                              hdf5_input_wall, &B_offload_array, NULL,
                              &plasma_offload_array, NULL, &wall_offload_array,
                              NULL, NULL);

    B_field_data B_data;
    B_field_init(&B_data, &sim.B_offload_data, B_offload_array);

    plasma_data plasma_data;
    plasma_init(&plasma_data, &sim.plasma_offload_data, plasma_offload_array);

    wall_data wall_data;
    wall_init(&wall_data, &sim.wall_offload_data, wall_offload_array);

    random_data rng;
    random_init(&rng, 12345);

    int n_inj;
    nbi_injector* inj;
    hid_t f = hdf5_open("ascot.h5");
    hdf5_nbi_read(f, &n_inj, &inj);

    for(int i=0; i < n_inj; i++) {
        printf("Injector %d:\n", i+1);
        printf("id: %d\n", inj[i].id);
        printf("n_beamlet: %d\n", inj[i].n_beamlet);
/*        printf("beamlet_x: ");
        for(int j=0; j < inj[i].n_beamlet; j++) {
            printf("%le ", inj[i].beamlet_x[j]);
        }
        printf("\n");*/
        printf("power: %le\n", inj[i].power);
        printf("energy: %le\n", inj[i].energy);
        printf("efrac: %le %le %le\n", inj[i].efrac[0], inj[i].efrac[1], inj[i].efrac[2]);
        printf("divergence: %le %le %le %le %le\n", inj[i].div_h, inj[i].div_v, inj[i].div_halo_frac, inj[i].div_halo_h, inj[i].div_halo_v);
        printf("anum: %d\n", inj[i].anum);
        printf("znum: %d\n", inj[i].znum);
        printf("\n");
    }

    int nprt = 100000;
    particle* p = (particle*) malloc(nprt*sizeof(particle));

    nbi_generate(nprt, p, &inj[0], &B_data, &plasma_data, &wall_data, &rng);

    FILE* of = fopen("nbi.out", "w");
    for(int i=0; i < nprt; i++) {
        fprintf(of,"%d %le %le %le %le %le %le\n", i, p[i].r, p[i].phi,
               p[i].z, p[i].v_r, p[i].v_phi, p[i].v_z);
    }
    fclose(of);

    return 0;
}

