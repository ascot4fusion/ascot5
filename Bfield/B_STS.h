/**
 * @file B_STS.h
 * @brief Header file for B_STS.c
 */
#ifndef B_STS_H
#define B_STS_H
#include "../ascot5.h"
#include "../symmetry.h"
#include "../linint/linint1D.h" /* for 1D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */

/**
 * @brief stellarator magnetic field parameters that will be offloaded to target
 */
typedef struct {
    int psigrid_n_r;             /**< number of r grid points in psi data */
    int psigrid_n_phi;           /**< number of phi grid points in psi data */
    int psigrid_n_z;             /**< number of z grid points in psi data */
    real psigrid_r_min;          /**< minimum r coordinate in the grid in psi data */
    real psigrid_r_max;          /**< maximum r coordinate in the grid in psi data */
    real psigrid_r_grid;         /**< r grid interval (r_max-r_min)/(n_r-1) in psi data */
    real psigrid_phi_min;        /**< minimum z coordinate in the grid in psi data */
    real psigrid_phi_max;        /**< maximum z coordinate in the grid in psi data */
    real psigrid_phi_grid;       /**< z grid interval (z_max-z_min)/(n_z-1) in psi data */
    real psigrid_z_min;          /**< minimum z coordinate in the grid in psi data */
    real psigrid_z_max;          /**< maximum z coordinate in the grid in psi data */
    real psigrid_z_grid;         /**< z grid interval (z_max-z_min)/(n_z-1) in psi data */
                                 
    int Bgrid_n_r;               /**< number of r grid points in B data */
    int Bgrid_n_phi;             /**< number of phi grid points in B data */
    int Bgrid_n_z;               /**< number of z grid points in B data */
    real Bgrid_r_min;            /**< minimum r coordinate in the grid in B data */
    real Bgrid_r_max;            /**< maximum r coordinate in the grid in B data */
    real Bgrid_r_grid;           /**< r grid interval (r_max-r_min)/(n_r-1) in B data */
    real Bgrid_phi_min;          /**< minimum phi coordinate in the grid in B data */
    real Bgrid_phi_max;          /**< maximum phi coordinate in the grid in B data */
    real Bgrid_phi_grid;         /**< phi grid interval 2pi/(n_phi-1) in B data */
    real Bgrid_z_min;            /**< minimum z coordinate in the grid in B data */
    real Bgrid_z_max;            /**< maximum z coordinate in the grid in B data */
    real Bgrid_z_grid;           /**< z grid interval (z_max-z_min)/(n_z-1) in B data */
                                 
    real period_length;          /**< length of one toroidal period in radians */
    symmetry_type symmetry_mode; /**< symmetry mode used */
                                 
    real psi0;                   /**< sqrt(psi) value at magnetic axis */
    real psi1;                   /**< sqrt(psi) value at separatrix */
                                 
    int n_axis;                  /**< number of phi grid points for magnetic axis */
    real axis_min;               /**< minimum phi coordinate in the magnetic axis grid */
    real axis_max;               /**< maximum phi coordinate in the magnetic axis grid */
    real axis_grid;              /**< phi grid interval 2pi/(n_phi-1) */
                                 
    int offload_array_length;    /**< number of elements in offload_array */
} B_STS_offload_data;

/**
 * @brief stellarator magnetic field parameters on the target
 */
typedef struct {
    real psi0;                   /**< sqrt(psi) value at magnetic axis */
    real psi1;                   /**< sqrt(psi) value at separatrix */
    real period_length;          /**< length of one toroidal period in radians */
    symmetry_type symmetry_mode; /**< symmetry mode used */
    linint1D_data axis_r;        /**< r coordinate of magnetic axis */
    linint1D_data axis_z;        /**< z coordinate of magnetic axis */
    interp3D_data psi;           /**< pointer to start of psi interpolation data struct */
    interp3D_data B_r;           /**< pointer to start of B_r interpolation data struct */
    interp3D_data B_phi;         /**< pointer to start of B_phi interpolation data struct */
    interp3D_data B_z;           /**< pointer to start of B_z interpolation data struct */
} B_STS_data;

int B_STS_init_offload(B_STS_offload_data* offload_data, real** offload_array);
void B_STS_free_offload(B_STS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
int B_STS_init(B_STS_data* Bdata, B_STS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_psi(real psi[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_rho(real rho[], real psi, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_B(real B[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_B_dB(real B_dB[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_get_axis_r(real* axis_r, B_STS_data* Bdata, real phi);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_get_axis_z(real* axis_z, B_STS_data* Bdata, real phi);
#pragma omp end declare target   
#endif
