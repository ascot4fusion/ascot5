/**
 * @file B_3DS.h
 * @brief Header file for B_3DS.c
 */
#ifndef B_3DS_H
#define B_3DS_H
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 3D magnetic field parameters on the host
 */
typedef struct {
    int psigrid_n_r;     /**< Number of R grid points in psi data             */
    int psigrid_n_z;     /**< Number of z grid points in psi data             */
    real psigrid_r_min;  /**< Minimum R grid point in psi data [m]            */
    real psigrid_r_max;  /**< Maximum R grid point in psi data [m]            */
    real psigrid_z_min;  /**< Minimum z grid point in psi data [m]            */
    real psigrid_z_max;  /**< Maximum z grid point in psi data [m]            */

    int Bgrid_n_r;       /**< Number of R grid points in B data               */
    int Bgrid_n_z;       /**< Number of z grid points in B data               */
    real Bgrid_r_min;    /**< Minimum R coordinate in the grid in B data [m]  */
    real Bgrid_r_max;    /**< Maximum R coordinate in the grid in B data [m]  */
    real Bgrid_z_min;    /**< Minimum z coordinate in the grid in B data [m]  */
    real Bgrid_z_max;    /**< Maximum z coordinate in the grid in B data [m]  */
    int Bgrid_n_phi;     /**< Number of phi grid points in B data             */
    real Bgrid_phi_min;  /**< Minimum phi grid point in B data [rad]          */
    real Bgrid_phi_max;  /**< Maximum phi grid point in B data [rad]          */

    real psi0;           /**< Poloidal flux value at magnetic axis [V*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    real axis_r;         /**< R coordinate of magnetic axis [m]               */
    real axis_z;         /**< z coordinate of magnetic axis [m]               */
    int offload_array_length; /**< Number of elements in offload_array        */
} B_3DS_offload_data;

/**
 * @brief 3D magnetic field parameters on the target
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [v*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    real axis_r;         /**< R coordinate of magnetic axis [m]               */
    real axis_z;         /**< z coordinate of magnetic axis [m]               */
    interp2D_data psi;   /**< 2D psi interpolation data struct                */
    interp3D_data B_r;   /**< 3D B_r interpolation data struct                */
    interp3D_data B_phi; /**< 3D B_phi interpolation data struct              */
    interp3D_data B_z;   /**< 3D B_z interpolation data struct                */
} B_3DS_data;

int B_3DS_init_offload(B_3DS_offload_data* offload_data, real** offload_array);
void B_3DS_free_offload(B_3DS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_3DS_init(B_3DS_data* Bdata, B_3DS_offload_data* offload_data,
                real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_psi(real psi[1], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_rho(real rho[1], real psi, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_B(real B[3], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DS_get_axis_r(B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DS_get_axis_z(B_3DS_data* Bdata);
#pragma omp end declare target
#endif
