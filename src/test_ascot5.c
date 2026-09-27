/**
 * @file test_ascot5.c
 * Contains the test program for running simulations with fixed inputs.
 */
#include "ascot.h"
#include "defines.h"
#include "datatypes.h"
#include <stdio.h>

/**
 * Run a simulation with fixed inputs.
 *
 * The purpose of this program is to aid in developing performance sensitive
 * parts of the code. Running this program directly instead of calling the
 * C library via Python allows for easier profiling and debugging of the C code.
 */
int main(void) {
    printf("test_ascot5\n");
    size_t nmrk = 1000;
    Simulation sim;
    State mrk[nmrk];
    ascot_solve_distribution(&sim, nmrk, mrk);
    return 0;
}