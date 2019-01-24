/**
 * @file B_3DST.h
 * @brief Header file for B_3DST.c
 */
#ifndef B_3DST_H
#define B_3DST_H

#include "../ascot5.h"
#include "../spline/interp2D.h" /* for 2D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */
//#include "../spline/interp4D.h"

/**
 * @brief 3D magnetic field parameters on the host
 */
typedef struct {
    real time[N_MAX_TIME_SLICE];

    int psigrid_n_r;            /**< number of r grid points in psi data */
    int psigrid_n_z;            /**< number of z grid points in psi data */
    int psigrid_n_phi;          /**< number of phi grid points in psi data */
    real psigrid_r_min;         /**< minimum r coordinate in the grid in psi data */
    real psigrid_r_max;         /**< maximum r coordinate in the grid in psi data */
    real psigrid_r_grid;        /**< r grid interval (r_max-r_min)/(n_r-1) in psi data */
    real psigrid_z_min;         /**< minimum z coordinate in the grid in psi data */
    real psigrid_z_max;         /**< maximum z coordinate in the grid in psi data */
    real psigrid_z_grid;        /**< z grid interval (z_max-z_min)/(n_z-1) in psi data */

    int Bgrid_n_r;              /**< number of r grid points in B data */
    int Bgrid_n_z;              /**< number of z grid points in B data */
    int n_phi;                  /**< number of phi grid points in B data */
    int n_time;                 /**< number of time grid points in B data */
    real Bgrid_r_min;           /**< minimum r coordinate in the grid in B data */
    real Bgrid_r_max;           /**< maximum r coordinate in the grid in B data */
    real Bgrid_r_grid;          /**< r grid interval (r_max-r_min)/(n_r-1) in B data */
    real Bgrid_z_min;           /**< minimum z coordinate in the grid in B data */
    real Bgrid_z_max;           /**< maximum z coordinate in the grid in B data */
    real Bgrid_z_grid;          /**< z grid interval (z_max-z_min)/(n_z-1) in B data */
    real phi_min;               /**< minimum phi coordinate in the grid in B data */
    real phi_max;               /**< maximum phi coordinate in the grid in B data */
    real phi_grid;              /**< phi grid interval 2pi/(n_phi-1) in B data */
    real time_min;              /**< minimum time coordinate in the grid in B data */
    real time_max;              /**< maximum time coordinate in the grid in B data */
    real time_grid;             /**< time grid interval (time_max-time_min)/(n_time-1) in B data */

    real psi0;                  /**< sqrt(psi) value at magnetic axis */
    real psi1;                  /**< sqrt(psi) value at separatrix */
    real axis_r;                /**< r coordinate of magnetic axis */
    real axis_z;                /**< z coordinate of magnetic axis */
    int offload_array_length;   /**< number of elements in offload_array */
} B_3DST_offload_data;

/**
 * @brief 3D magnetic field parameters on the target
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [v*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    real axis_r;         /**< R coordinate of magnetic axis [m]               */
    real axis_z;         /**< z coordinate of magnetic axis [m]               */
    interp2D_data psi;   /**< 2D psi interpolation data struct                */
    interp4D_data B_r;   /**< 3D B_r interpolation data struct                */
    interp4D_data B_phi; /**< 3D B_phi interpolation data struct              */
    interp4D_data B_z;   /**< 3D B_z interpolation data struct                */
} B_3DST_data;

void B_3DST_init_offload(B_3DST_offload_data* offload_data, real** offload_array);
void B_3DST_free_offload(B_3DST_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_3DST_init(B_3DST_data* Bdata, B_3DST_offload_data* offload_data,
                real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_3DST_eval_psi(real psi[1], real r, real phi, real z, B_3DST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DST_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
			   B_3DST_data* Bdata); //this is not used really
#pragma omp declare simd uniform(Bdata)
a5err B_3DST_eval_rho(real rho[1], real psi, B_3DST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DST_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_3DST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DST_eval_B(real B[3], real r, real phi, real z, real time, B_3DST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DST_eval_B_dB(real B_dB[12], real r, real phi, real z, real time,
                      B_3DST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DST_get_axis_r(B_3DST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DST_get_axis_z(B_3DST_data* Bdata);
#pragma omp end declare target
#endif
