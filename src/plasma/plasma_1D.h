/**
 * @file plasma_1D.h
 * @brief Header file for plasma_1D.c
 */
#ifndef PLASMA_1D_H
#define PLASMA_1D_H
#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct {
    int n_rho;     /**< number of rho values in the data             */
    int n_species; /**< number of plasma species including electrons */
    real* mass;    /**< plasma species masses [kg]                   */
    real* charge;  /**< plasma species charges [C]                   */
    int* anum;     /**< ion species atomic number                    */
    int* znum;     /**< ion species charge number                    */
    real* rho;     /**< pointer to start of rho values               */
    real* temp;    /**< pointer to start of temperatures             */
    real* dens;    /**< pointer to start of densities                */
    real* vtor;    /**< pointer to start of toroidal rotation        */
} plasma_1D_data;

int plasma_1D_init(plasma_1D_data* data, int nrho, int nion, real* rho,
                   int* anum, int* znum, real* mass, real* charge,
                   real* Te, real* Ti, real* ne, real* ni);
void plasma_1D_free(plasma_1D_data* data);
void plasma_1D_offload(plasma_1D_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1D_eval_temp(real* dens, real rho, int species,
                          plasma_1D_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1D_eval_dens(real* temp, real rho, int species,
                          plasma_1D_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1D_eval_densandtemp(real* dens, real* temp, real rho,
                                 plasma_1D_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1D_eval_flow(real* vflow, real rho, plasma_1D_data* pls_data);
DECLARE_TARGET_END

#endif
