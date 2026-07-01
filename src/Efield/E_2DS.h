/**
 * @file E_2DS.h
 * @brief Header file for E_2DS.c
 *
 * Contains declaration of E_2DS_field_data struct.
 */
#ifndef E_2DS_H
#define E_2DS_H
#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 2D spline electric field parameters
 */
typedef struct {
    interp2D_data vpot;  /**< Interpolation struct for the 2D potential */
} E_2DS_data;

int E_2DS_init(E_2DS_data* data, int nr, real rmin, real rmax, int nz,
    real zmin, real zmax, real* vpot);
void E_2DS_free(E_2DS_data* data);
void E_2DS_offload(E_2DS_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Edata)
a5err E_2DS_eval_E(real E[3], real r, real z, E_2DS_data* Edata);
DECLARE_TARGET_END
#endif