/**
 * @file B_3DN.h
 * @brief Header file for B_3DN.c
 */
#ifndef B_3DN_H
#define B_3DN_H
#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"     // TODO: REMOVE WHEN DONE
#include "../neural_network/neural_network.h"

/**
 * @brief 3D magnetic field parameters
 */

// TODO: REPLACE WITH THE WEIGHTS
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [v*s*m^-1] */   //THIS WILL STAY
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */   //THIS WILL STAY
    real axis_r;         /**< R coordinate of magnetic axis [m]               */   //THIS WILL STAY
    real axis_z;         /**< z coordinate of magnetic axis [m]               */   //THIS WILL STAY
    neural2D_data psi_neural_data;   /**< 2D psi interpolation data struct                */   // OK ?
    neural3D_data B_rphiz_neural_data;   /**< 3D B_r interpolation data struct                */   // OK ?
    //neural3D_data B_phi; /**< 3D B_phi interpolation data struct              */   // OK ?
    //neural3D_data B_z;   /**< 3D B_z interpolation data struct                */   // OK ?
} B_3DN_data;

int B_3DN_init(B_3DN_data* data,
               int p_n_r, real p_r_min, real p_r_max,
               int p_n_z, real p_z_min, real p_z_max,
               int b_n_r, real b_r_min, real b_r_max,
               int b_n_phi, real b_phi_min, real b_phi_max,
               int b_n_z, real b_z_min, real b_z_max,
               real axis_r, real axis_z, real psi0, real psi1,
               real* psi, real* B_r, real* B_phi, real* B_z);
void B_3DN_free(B_3DN_data* data);
void B_3DN_offload(B_3DN_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_3DN_eval_psi(real* psi, real r, real phi, real z, B_3DN_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_3DN_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_3DN_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_3DN_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_3DN_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_3DN_eval_B(real B[3], real r, real phi, real z, B_3DN_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_3DN_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_3DN_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_3DN_get_axis_rz(real rz[2], B_3DN_data* Bdata);
DECLARE_TARGET_END
#endif
