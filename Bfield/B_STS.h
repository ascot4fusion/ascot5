/**
 * @file B_STS.h
 * @brief Header file for B_STS.c
 */
#ifndef B_STS_H
#define B_STS_H
#include "../ascot5.h"
#include "../linint/linint.h"
#include "../spline/interp.h"

/**
 * @brief stellarator magnetic field parameters on the host
 */
typedef struct {
    int psigrid_n_r;       /**< Number of R grid points in psi data           */
    int psigrid_n_z;       /**< Number of z grid points in psi data           */
    int psigrid_n_phi;     /**< Number of phi grid points in psi data         */
    real psigrid_r_min;    /**< Minimum R grid point in psi data [m]          */
    real psigrid_r_max;    /**< Maximum R grid point in psi data [m]          */
    real psigrid_z_min;    /**< Minimum z grid point in psi data [m]          */
    real psigrid_z_max;    /**< Maximum z grid point in psi data [m]          */
    real psigrid_phi_min;  /**< Minimum phi grid point in psi data [rad]      */
    real psigrid_phi_max;  /**< Maximum phi grid point in psi data [rad]      */

    int Bgrid_n_r;       /**< Number of R grid points in B data               */
    int Bgrid_n_z;       /**< Number of z grid points in B data               */
    int Bgrid_n_phi;     /**< Number of phi grid points in B data             */
    real Bgrid_r_min;    /**< Minimum R coordinate in the grid in B data [m]  */
    real Bgrid_r_max;    /**< Maximum R coordinate in the grid in B data [m]  */
    real Bgrid_z_min;    /**< Minimum z coordinate in the grid in B data [m]  */
    real Bgrid_z_max;    /**< Maximum z coordinate in the grid in B data [m]  */
    real Bgrid_phi_min;  /**< Minimum phi grid point in B data [rad]          */
    real Bgrid_phi_max;  /**< Maximum phi grid point in B data [rad]          */

    real psi0;           /**< Poloidal flux value at magnetic axis [V*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    int offload_array_length; /**< Number of elements in offload_array        */

    int n_axis;          /**< Number of phi grid points in axis data          */
    real axis_min;       /**< Minimum phi grid point in axis data [rad]       */
    real axis_max;       /**< Maximum phi grid point in axis data [rad]       */
    real axis_grid;      /**< phi grid interval in axis data [rad]            */
} B_STS_offload_data;

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

int B_STS_init_offload(B_STS_offload_data* offload_data, real** offload_array);
void B_STS_free_offload(B_STS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_STS_init(B_STS_data* Bdata, B_STS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_psi(real* psi, real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_rho(real* rho, real psi, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_B(real B[3], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_get_axis_r(real* axis_r, B_STS_data* Bdata, real phi);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_get_axis_z(real* axis_z, B_STS_data* Bdata, real phi);
#pragma omp end declare target
#endif
