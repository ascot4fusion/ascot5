/**
 * @file plasma_1Dt.h
 * @brief Header file for plasma_1Dt.c
 */
#ifndef PLASMA_1DT_H
#define PLASMA_1DT_H
#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct {
    int n_rho;     /**< number of rho values in the data                */
    int n_time;    /**< number of time points                           */
    int n_species; /**< number of plasma species including electrons    */
    real* mass;    /**< plasma species masses [kg]                      */
    real* charge;  /**< plasma species charges [C]                      */
    int* anum;     /**< ion species atomic number                       */
    int* znum;     /**< ion species charge number                       */
    real* rho;     /**< pointer to start of rho values                  */
    real* time;    /**< pointer to start of time values                 */
    real* temp;    /**< pointer to start of temperatures                */
    real* dens;    /**< pointer to start of densities                   */
} plasma_1Dt_data;

int plasma_1Dt_init(plasma_1Dt_data* data, int nrho, int ntime, int nion,
                    real* rho, real* time, int* anum, int* znum, real* mass,
                    real* charge, real* Te, real* Ti, real* ne, real* ni);
void plasma_1Dt_free(plasma_1Dt_data* pls_data);
void plasma_1Dt_offload(plasma_1Dt_data* pls_data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1Dt_eval_temp(real* dens, real rho, real t, int species,
                          plasma_1Dt_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1Dt_eval_dens(real* temp, real rho, real t, int species,
                          plasma_1Dt_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_1Dt_eval_densandtemp(real* dens, real* temp, real rho, real t,
                                 plasma_1Dt_data* pls_data);
DECLARE_TARGET_END

#endif
