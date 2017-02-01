/**
 * @file ascotpy.c
 * @brief Python interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "B_field.h"

static real* B_offload_array;
static B_field_data B_data;

void ascotpy_isPointInsideWall(int handle, double* locRPZ, int* isInside,
                               int* err) {
}

void ascotpy_initialize(char* folder, int handle, int* err) {
    B_field_offload_data B_offload_data;
    B_field_init_offload(&B_offload_data, &B_offload_array);
    B_field_init(&B_data, &B_offload_data, B_offload_array);
}

void ascotpy_free(int handle, int* err) {

}

void ascotpy_evalB(int handle, double* locRPZ, double* BRPZ, int* err) {
    real B[3][NSIMD];
    B_field_eval_B(0, B, locRPZ[0], locRPZ[1], locRPZ[2], &B_data);
    BRPZ[0] = B[0][0];
    BRPZ[1] = B[1][0];
    BRPZ[2] = B[2][0];
}

void ascotpy_evalBjacB(int handle, double* locRPZ, double* BRPZ, double* jacB,
                  int* err) {
    real B_dB[12][NSIMD];
    B_field_eval_B_dB(0, B_dB, locRPZ[0], locRPZ[1], locRPZ[2], &B_data);
    BRPZ[0] = B_dB[0][0];
    BRPZ[1] = B_dB[4][0];
    BRPZ[2] = B_dB[8][0];
    jacB[0] = B_dB[1][0];
    jacB[1] = B_dB[2][0];
    jacB[2] = B_dB[3][0];
    jacB[3] = B_dB[5][0];
    jacB[4] = B_dB[6][0];
    jacB[5] = B_dB[7][0];
    jacB[6] = B_dB[9][0];
    jacB[7] = B_dB[10][0];
    jacB[8] = B_dB[11][0];
}
