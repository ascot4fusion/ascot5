/**
 * @file B_STS.h
 * @brief Header file for B_STS.c
 */
#ifndef B_STS_H
#define B_STS_H
#include "../offload.h"
#include "../ascot5.h"
#include "../linint/linint.h"
#include "../spline/interp.h"

/**
 * @brief stellarator magnetic field parameters on the target
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [v*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    linint1D_data axis_r;/**< 1D axis r-value interpolation data struct       */
    linint1D_data axis_z;/**< 1D axis z-value interpolation data struct       */
    interp3D_data psi;   /**< 3D psi interpolation data struct                */
    interp3D_data B_r;   /**< 3D B_r interpolation data struct                */
    interp3D_data B_phi ;/**< 3D B_phi interpolation data struct              */
    interp3D_data B_z;   /**< 3D B_z interpolation data struct                */
} B_STS_data;

int B_STS_init(B_STS_data* data,
               int p_n_r, real p_r_min, real p_r_max,
               int p_n_phi, real p_phi_min, real p_phi_max,
               int p_n_z, real p_z_min, real p_z_max,
               int b_n_r, real b_r_min, real b_r_max,
               int b_n_phi, real b_phi_min, real b_phi_max,
               int b_n_z, real b_z_min, real b_z_max,
               int naxis, real axis_min, real axis_max,
               real* axis_r, real* axis_z, real psi0, real psi1,
               real* psi, real* B_r, real* B_phi, real* B_z);
void B_STS_free(B_STS_data* data);
void B_STS_offload(B_STS_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_STS_eval_psi(real* psi, real r, real phi, real z, B_STS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_STS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_STS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_STS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_STS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_STS_eval_B(real B[3], real r, real phi, real z, B_STS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_STS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_STS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_STS_get_axis_rz(real rz[2], B_STS_data* Bdata, real phi);
DECLARE_TARGET_END
#endif
