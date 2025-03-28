/**
 * @file B_GS.h
 * @brief Header file for B_GS.c
 */
#ifndef B_GS_H
#define B_GS_H
#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"

/**
 * @brief Analytic magnetic field parameters on the target
 */
typedef struct {
    real R0;                  /**< Major radius R coordinate                  */
    real z0;                  /**< Midplane z coordinate                      */
    real raxis;               /**< Magnetic axis R coordinate                 */
    real zaxis;               /**< Magnetic axis z coordinate                 */
    real B_phi0;              /**< On-axis toroidal field                     */
    real psi0;                /**< Poloidal flux at axis [V*s*m^-1]           */
    real psi1;                /**< Poloidal flux at separatrix [V*s*m^-1]     */
    real psi_mult;            /**< Psi multiplier                             */
    real psi_coeff[14];       /**< Coefficients for evaluating psi
                                   [c_1, c_2, ..., c_13, A]                   */
    int Nripple;              /**< Number of toroidal field coils             */
    real a0;                  /**< Minor radius [m]                           */
    real alpha0;              /**< Ripple r-dependency, delta ~ (r/a0)^alpha0 */
    real delta0;              /**< Ripple strength                            */
} B_GS_data;

int B_GS_init(B_GS_data* data, real R0, real z0, real raxis, real zaxis,
              real B_phi0, real psi0, real psi1, real psi_mult, real c[14],
              int Nripple, real a0, real alpha0, real delta0);
void B_GS_free(B_GS_data* data);
void B_GS_offload(B_GS_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_GS_eval_B(real B[3], real r, real phi, real z, B_GS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_GS_eval_psi(real* psi, real r, real phi, real z, B_GS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_GS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                         B_GS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_GS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                         B_GS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_GS_eval_B_dB(real B_dB[12], real r, real phi, real z, B_GS_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_GS_get_axis_rz(real rz[2], B_GS_data* Bdata);
DECLARE_TARGET_END

#endif
