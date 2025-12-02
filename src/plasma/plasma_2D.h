/**
 * @file plasma_2D.h
 * @brief Header file for plasma_2D.c
 */
#ifndef PLASMA_2D_H
#define PLASMA_2D_H
#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"
#include "../linint/linint.h"

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct {
    int n_species;       /**< number of plasma species including electrons */
    real* mass;          /**< plasma species masses [kg]                   */
    real* charge;        /**< plasma species charges [C]                   */
    int* anum;           /**< ion species atomic number                    */
    int* znum;           /**< ion species charge number                    */
    linint2D_data* temp; /**< temperature interpolant                      */
    linint2D_data* dens; /**< density interpolant                          */
    linint2D_data* vtor; /**< toroidal rotation interpolant                */
} plasma_2D_data;

int plasma_2D_init(plasma_2D_data* data, int nr, int nz, int nion,
                   real r_min, real r_max, real z_min, real z_max,
                   int* anum, int* znum, real* mass, real* charge,
                   real* Te, real* Ti, real* ne, real* ni, real* vtor);
void plasma_2D_free(plasma_2D_data* data);
void plasma_2D_offload(plasma_2D_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_2D_eval_temp(real* dens, real r, real z, int species,
                          plasma_2D_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_2D_eval_dens(real* temp, real r, real z, int species,
                          plasma_2D_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_2D_eval_densandtemp(real* dens, real* temp, real r, real z,
                                 plasma_2D_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_2D_eval_flow(real* vflow, real r, real z,
                          plasma_2D_data* pls_data);
DECLARE_TARGET_END

#endif