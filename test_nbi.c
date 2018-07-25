/**
 * @file test_nbi.c
 * @brief Test program for NBI functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "nbi.h"
#include "hdf5io/hdf5_nbi.h"
#include "hdf5io/hdf5_helpers.h"
#include "suzuki.h"

int main(int argc, char** argv) {
/*    int Z[] = {1,4};
    int A[] = {3,9};
    real ni[] = {5e19, 1e17};
    printf("%le\n", suzuki_sigmav(150, 5.0025e19, 5000, 2, ni, A, Z));

    return 0;*/

    hdf5_init();
    hid_t f = hdf5_open("ascot.h5");

    int n_inj;
    nbi_injector* inj;
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
        printf("divergence: %le %le %le\n", inj[i].divergence[0], inj[i].divergence[1], inj[i].divergence[2]);
        printf("anum: %d\n", inj[i].anum);
        printf("znum: %d\n", inj[i].znum);
        printf("\n");
    }

    return 0;
}

