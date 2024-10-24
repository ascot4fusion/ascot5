/**
 * @file B_2DS.h
 * @brief Header file for B_2DS.c
 *
 * Contains declaration of B_2DS_offload_data and B_2DS_data structs.
 */
#ifndef B_2DS_H
#define B_2DS_H
#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 2D magnetic field parameters
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

int B_2DS_init(B_2DS_data* data,
               int n_r, real r_min, real r_max,
               int n_z, real z_min, real z_max,
               real axis_r, real axis_z, real psi0, real psi1,
               real* psi, real* B_r, real* B_phi, real* B_z);
void B_2DS_free(B_2DS_data* data);
void B_2DS_offload(B_2DS_data* data);
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
