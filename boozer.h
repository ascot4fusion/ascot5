/**
 * @file boozer.h
 * @brief Header file for boozer.c
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "ascot5.h"
#include "error.h"
#include "spline/interp.h"

/**
 * @brief offload data for coordinate maps for moving from (R,phi,z) to 
 *        boozer (psi,theta,zeta)
 *
 */
typedef struct {
    int  nr;       /**< Number of R grid points for psi_rz                    */
    real rmin;     /**< Minimum R for psi_rz                                  */
    real rmax;     /**< Maximum R for psi_rz                                  */
    int  nz;       /**< Number of z grid points for psi_rz                    */
    real zmin;     /**< Minimum z for psi_rz                                  */
    real zmax;     /**< Maximum z for psi_rz                                  */
    int  npsi;     /**< Number of psi grid points for other fields            */
    real psimin;   /**< Minimum psi in other fields                           */
    real psimax;   /**< Maximum psi in other fields                           */
    real psi_inner;/**< psi at the inner edge (center of plasma) of the grid  */
    real psi_outer;/**< psi at the outer edge (~separatrix) of the grid       */
    int  ntheta;   /**< number of theta angles (both boozer and geometric)    */
    real thetamin; /**< minimum theta (both the boozer and geometric)         */
    real thetamax; /**< maximum theta (both the boozer and geometric)         */
    real r0;       /**< R location of the axis for defining geometric theta   */
    real z0;       /**< z location of the axis for defining geometric theta   */
    int  nrzs;     /**< number of elements in rs and zs                       */
    int  offload_array_length; /**< Number of elements in offload_array       */
} boozer_offload_data;

/**
 * @brief Boozer parameters on the target
 */
typedef struct {
    real rmin;     /**< Minimum R for psi_rz */
    real rmax;     /**< Maximum R for psi_rz */
    real zmin;     /**< Minimum z for psi_rz */
    real zmax;     /**< Maximum z for psi_rz */
    real psimin;   /**< Minimum psi in other fields */
    real psimax;   /**< Maximum psi in other fields */
    real psi_inner;/**< psi at the inner edge (center of plasma) of the grid  */
    real psi_outer;/**< psi at the outer edge (~separatrix) of the grid       */
    real thetamin; /**< minimum theta (both the boozer and geometric) */
    real thetamax; /**< maximum theta (both the boozer and geometric) */
    real r0;       /**< R location of the axis for defining geometric theta */
    real z0;       /**< z location of the axis for defining geometric theta */
    real* rs;      /**< R points of outermost poloidal psi-surface contour,
		        nrzs elements, the first and last points are the same */
    real* zs;      /**< z points of outermost poloidal psi-surface contour,
		        nrzs elements, the first and last points are the same */
    int  nrzs;     /**< number of elements in rs and zs */
    interp2D_data nu_psitheta; /**< the nu-function, phi=zeta+nu(psi,theta),
                                    with phi the cylindrical angle            */
    interp2D_data theta_psithetageom; /**< boozer_theta(psi,theta_geometric)  */
    interp2D_data psi_rz; /**< psi(R,z) */
} boozer_data;

int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array);
void boozer_free_offload(boozer_offload_data* offload_data,
                         real** offload_array);

#pragma omp declare target
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array);

#pragma omp declare simd uniform(boozerdata)
a5err boozer_eval_psithetazeta(real psithetazeta[12], real r, real phi, real z,
			       boozer_data* boozerdata);

#pragma omp declare simd uniform(boozerdata)
a5err boozer_eval_psinormalized(real psi[4], real psin[4],
				boozer_data* boozerdata);

#pragma omp end declare target

#endif
