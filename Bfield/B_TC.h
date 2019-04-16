/**
 * @file B_TC.h
 * @brief Header file for B_TC.c
 *
 * Contains declaration of B_TC_offload_data and B_TC_data structs.
 */
#ifndef B_TC_H
#define B_TC_H
#include "../ascot5.h"
#include "../error.h"

/**
 * @brief TC magnetic field parameters that will be offloaded to target
 */
typedef struct {
    real axisr;  /**< A value that is returned when magnetic axis
                      R coordinate is requested [m]                           */
    real axisz;  /**< A value that is returned when magnetic axis
                      z coordinate is requested [m]                           */
    real psival; /**< A value that is returned when poloidal magnetic flux
                      value is requested [V*s*m^-1]                           */
    real rhoval; /**< A value that is returened when normalized poloidal flux
                      value is requested                                      */
    real B[3];   /**< Magnetic field at origo: [B_x, B_y, B_z] [T]            */
    real dB[9];  /**< Magnetic field Jacobian:
                      [dB_x/dx, dB_x/dy, dB_x/dz,
                       dB_y/dx, dB_y/dy, dB_y/dz,
                       dB_z/dx, dB_z/dy, dB_z/dz] [T/m]                       */
    int offload_array_length; /**< Number of elements in offload_array; zero
                                   always for this magnetic field type        */
} B_TC_offload_data;

/**
 * @brief TC magnetic field parameters on the target
 */
typedef struct {
    real axisr;  /**< A value that is returned when magnetic axis
                      R coordinate is requested [m]                           */
    real axisz;  /**< A value that is returned when magnetic axis
                      z coordinate is requested [m]                           */
    real psival; /**< A value that is returned when poloidal magnetic flux
                      value is requested [V*s*m^-1]                           */
    real rhoval; /**< A value that is returened when normalized poloidal flux
                      value is requested                                      */
    real* B;     /**< Magnetic field at origo: [B_x, B_y, B_z] [T]            */
    real* dB;    /**< Magnetic field Jacobian:
                      [dB_x/dx, dB_x/dy, dB_x/dz,
                       dB_y/dx, dB_y/dy, dB_y/dz,
                       dB_z/dx, dB_z/dy, dB_z/dz] [T/m]                       */
} B_TC_data;

int B_TC_init_offload(B_TC_offload_data* offload_data, real** offload_array);
void B_TC_free_offload(B_TC_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_TC_init(B_TC_data* Bdata, B_TC_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_TC_eval_B(real B[3], real r, real phi, real z, B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_TC_eval_psi(real* psi, real r, real phi, real z, B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_TC_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                         B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_TC_eval_rho(real* rho, real psi, B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_TC_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                         B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_TC_eval_B_dB(real B_dB[12], real r, real phi, real z, B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_TC_get_axis_r(B_TC_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_TC_get_axis_z(B_TC_data* Bdata);
#pragma omp end declare target

#endif
