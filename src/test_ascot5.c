/**
 * @file test_ascot5.c
 * Contains the test program for running simulations with fixed inputs.
 */
#include "ascot.h"
#include "defines.h"
#include "datatypes.h"
#include <stdio.h>

#include <stdint.h>
#include <philox.h>

/**
 * Run a simulation with fixed inputs.
 *
 * The purpose of this program is to aid in developing performance sensitive
 * parts of the code. Running this program directly instead of calling the
 * C library via Python allows for easier profiling and debugging of the C code.
 */
int main(void) {
    printf("test_ascot5\n");
    //size_t nmrk = 1000;
    //Simulation sim;
    //State mrk[nmrk];
    //ascot_solve_distribution(&sim, nmrk, mrk);
    philox4x32_ctr_t counter = {{0, 0, 0, 0}};
    philox4x32_ukey_t user_key = {{12345, 0}};
    philox4x32_key_t key = philox4x32keyinit(user_key);

    /* Generate four 32-bit pseudorandom integers. */
    philox4x32_ctr_t result = philox4x32(counter, key);

    for (int i = 0; i < 4; ++i) {
        printf("%u\n", (unsigned)result.v[i]);
    }
    return 0;
}
