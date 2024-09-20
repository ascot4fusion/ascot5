/**
 * @file plasma_1DS.h
 * @brief Header file for plasma_1DS.c
 */
#ifndef PLASMA_1DS_H
#define PLASMA_1DS_H
#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 1D spline plasma parameters on the target
 */
typedef struct {
    int n_species; /**< number of plasma species including electrons */
    real* mass;    /**< plasma species masses (kg)                   */
    real* charge;  /**< plasma species charges (C)                   */
    int* anum;     /**< ion species atomic number                    */
    int* znum;     /**< ion species charge number                    */
    interp1D_data temp[2]; /**< electron and ion temperature interpolation    */
    interp1D_data* dens;   /**< electron and ion density interpolation structs*/
} plasma_1DS_data;

int plasma_1DS_init(plasma_1DS_data* data, int nrho, real rhomin, real rhomax,
                    int nion, int* anum, int* znum, real* mass, real* charge,
                    real* Te, real* Ti, real* ne, real* ni);
void plasma_1DS_free(plasma_1DS_data* data);
void plasma_1DS_offload(plasma_1DS_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1DS_eval_temp(real* temp, real rho, int species,
                           plasma_1DS_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1DS_eval_dens(real* dens, real rho, int species,
                           plasma_1DS_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1DS_eval_densandtemp(real* dens, real* temp, real rho,
                                  plasma_1DS_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1DS_eval_rotation(real* vr, real* vphi, real* vz, real rho, real r,
                               plasma_1DS_data* pls_data);
DECLARE_TARGET_END

#endif
