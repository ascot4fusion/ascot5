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
 * @brief 3D magnetic field parameters when using neural networks.
 *
 * Psi comes from a (R,z) --> psi NN. B_rphiz comes from a (R, phi, z) --> (Br, Bphi, Bz) NN.
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [v*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    real axis_r;         /**< R coordinate of magnetic axis [m]               */
    real axis_z;         /**< z coordinate of magnetic axis [m]               */
    neural2D_data psi_neural_data;   /**< 2D (R,z) --> psi interpolation data struct                */
    neural3D_data B_rphiz_neural_data;   /**< 3D (R, phi, z) --> (Br, Bphi, Bz) interpolation data struct                */
} B_3DN_data;

int B_3DN_init(B_3DN_data* data,
               real axis_r, real axis_z, real psi0, real psi1,
               real* psi_weights, real* B_rphiz_weights, int psi_n_layers,
               int* psi_layer_dimensions, int B_rphiz_n_layers,
               int* B_rphiz_layer_dimensions,
               real psi_r_mean, real psi_z_mean, real psi_psi_mean,
               real psi_r_std, real psi_z_std, real psi_psi_std,
               real bfield_r_mean, real bfield_phi_mean, real bfield_z_mean,
               real bfield_Br_mean, real bfield_Bphi_mean, real bfield_Bz_mean,
               real bfield_r_std, real bfield_phi_std, real bfield_z_std,
               real bfield_Br_std, real bfield_Bphi_std, real bfield_Bz_std);
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
