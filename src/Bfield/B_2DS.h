/**
 * @file B_2DS.h
 * @brief Header file for B_2DS.c
 *
 * Contains declaration of B_2DS_offload_data and B_2DS_data structs.
 */
#ifndef B_2DS_H
#define B_2DS_H
#include "../offload_acc_omp.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 2D magnetic field parameters that will be offloaded to target
 */
typedef struct {
    int n_r;                  /**< Number of r grid points                    */
    int n_z;                  /**< Number of z grid points                    */
    real r_min;               /**< Minimum R coordinate in the grid [m]       */
    real r_max;               /**< Maximum R coordinate in the grid [m]       */
    real z_min;               /**< Minimum z coordinate in the grid [m]       */
    real z_max;               /**< Maximum z coordinate in the grid [m]       */
    real psi0;                /**< Poloidal flux at magnetic axis [V*s*m^-1]  */
    real psi1;                /**< Poloidal flux at separatrix [V*s*m^-1]     */
    real axis_r;              /**< R coordinate of magnetic axis [m]          */
    real axis_z;              /**< z coordinate of magnetic axis [m]          */
    int offload_array_length; /**< Number of elements in offload_array        */
} B_2DS_offload_data;

/**
 * @brief 2D magnetic field parameters on the target
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [V*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    real axis_r;         /**< R coordinate of magnetic axis [m]               */
    real axis_z;         /**< z coordinate of magnetic axis [m]               */
    interp2D_data psi;   /**< psi interpolation 2D spline struct              */
    interp2D_data B_r;   /**< B_R interpolation 2D spline struct              */
    interp2D_data B_phi; /**< B_phi interpolation 2D spline struct            */
    interp2D_data B_z;   /**< B_z interpolation 2D spline struct              */
} B_2DS_data;

int B_2DS_init_offload(B_2DS_offload_data* offload_data, real** offload_array);
void B_2DS_free_offload(B_2DS_offload_data* offload_data, real** offload_array);

void B_2DS_init(B_2DS_data* Bdata, B_2DS_offload_data* offload_data,
                real* offload_array);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_2DS_eval_psi(real* psi, real r, real phi, real z, B_2DS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_2DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_2DS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_2DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_2DS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_2DS_eval_B(real B[3], real r, real phi, real z, B_2DS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_2DS_eval_B_dB(real B_dB[12], real r, real phi, real z, B_2DS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_2DS_get_axis_rz(real rz[2], B_2DS_data* Bdata);
DECLARE_TARGET_END
#endif
