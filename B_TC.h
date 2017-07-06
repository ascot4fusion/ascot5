/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file B_TC.h
 * @brief Header file for B_TC.c
 */
#ifndef B_TC_H
#define B_TC_H
#include "ascot5.h"

/**
 * @brief Trivial Cartesian magnetic field parameters that will be offloaded to target
 */
typedef struct {
    real axisr;
    real axisz;
    real psival;
    real rhoval;
    int offload_array_length; /**< number of elements in offload_array */
} B_TC_offload_data;

/**
 * @brief Trivial Cartesian magnetic field parameters on the target
 */
typedef struct {
    real axisr;
    real axisz;
    real psival;
    real rhoval;
    real* B;                   /**< Magnetic field at origo */
    real* dB;                  /**< Magnetic field gradient */
} B_TC_data;

void B_TC_init_offload(B_TC_offload_data* offload_data, real** offload_array);
void B_TC_free_offload(B_TC_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_TC_init(B_TC_data* Bdata, B_TC_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
void B_TC_eval_B(real* B, real r, real phi, real z, B_TC_data* Bdata);
void B_TC_eval_psi(real* psi, real r, real phi, real z, B_TC_data* Bdata);
void B_TC_eval_psi_dpsi(real* psi, real r, real phi, real z, B_TC_data* Bdata);
void B_TC_eval_rho(real* rho, real psi, B_TC_data* Bdata);
void B_TC_eval_rho_drho(real* rho, real r, real phi, real z, B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_TC_eval_B_dB(real* B_dB, real r, real phi, real z, B_TC_data* Bdata);
real B_TC_get_axis_r(B_TC_data* Bdata);
real B_TC_get_axis_z(B_TC_data* Bdata);
#pragma omp end declare target

#endif
