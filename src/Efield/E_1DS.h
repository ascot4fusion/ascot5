/**
 * @author Joona Kontula joona.kontula@aalto.fi
 * @file E_1DS.h
 * @brief Header file for E_1DS.c
 *
 * Contains declaration of E_1DS_field_offload_data and E_1DS_field_data
 * structs.
 */
#ifndef E_1DS_H
#define E_1DS_H
#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"
#include "../B_field.h"

/**
 * @brief 1D spline electric field parameters on the target
 */
typedef struct {
    interp1D_data dV;  /**< dV_drho 1D linear interpolation struct */
} E_1DS_data;

int E_1DS_init(E_1DS_data* data, int nrho, real rhomin, real rhomax, real reff,
               real* dvdrho);
void E_1DS_free(E_1DS_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Edata,Bdata)
a5err E_1DS_eval_E(real E[3], real r, real phi, real z, E_1DS_data* Edata,
                   B_field_data* Bdata);
DECLARE_TARGET_END
#endif
