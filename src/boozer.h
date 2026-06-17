/**
 * @file boozer.h
 * @brief Header file for boozer.c
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "spline/interp.h"

/**
 * @brief Data for mapping between the cylindrical and Boozer coordinates
 */
typedef struct {
    real psi_min;  /**< Minimum psi in other fields                           */
    real psi_max;  /**< Maximum psi in other fields                           */
    real* rs;  /**< R points of outermost poloidal psi-surface contour,
                    nrzs elements, the first and last points are the same     */
    real* zs;  /**< z points of outermost poloidal psi-surface contour,
                    nrzs elements, the first and last points are the same     */
    int  nrzs; /**< number of elements in rs and zs                           */
    interp2D_data nu_psitheta; /**< the nu-function, phi=zeta+nu(psi,theta),
                                    with phi the cylindrical angle            */
    interp2D_data theta_psithetageom; /**< boozer_theta(psi,thetag)           */
} boozer_data;

int boozer_init(boozer_data* data, int npsi, real psi_min, real psi_max,
                int ntheta, int nthetag, real* nu, real* theta,
                int nrzs, real* rs, real* zs);
void boozer_free(boozer_data* data);
void boozer_offload(boozer_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata, boozerdata)
a5err boozer_eval_psithetazeta(real psithetazeta[12], int* isinside, real r,
                               real phi, real z, B_field_data* Bdata,
                               boozer_data* boozerdata);


#endif
