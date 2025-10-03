/**
 * @file boozer.h
 * Module for transforming between cylindrical and Boozer coordinates.
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "B_field.h"
#include "ascot5.h"
#include "error.h"
#include "interp.h"

/**
 * Data for mapping between the cylindrical and Boozer coordinates.
 */
typedef struct
{
    real psi_min; /**< Minimum psi in other fields.                           */
    real psi_max; /**< Maximum psi in other fields.                           */
    real *rs;     /**< R points of outermost poloidal psi-surface contour,
                       nrzs elements, the first and last points are the same.     */
    real *zs;     /**< z points of outermost poloidal psi-surface contour,
                       nrzs elements, the first and last points are the same.     */
    int nrzs; /**< Number of elements in rs and zs.                           */
    interp2D_data nu_psitheta; /**< The nu-function, phi=zeta+nu(psi,theta),
                                    with phi the cylindrical angle. */
    interp2D_data theta_psithetageom; /**< Boozer_theta(psi,thetag). */
} boozer_data;

/**
 * Initialize boozer coordinate transformation.
 *
 * Multidimensional arrays must be stored as
 * - nu(psi_i, theta_j)     = array[j*npsi + i]
 * - theta(psi_i, thetag_j) = array[j*npsi + i]
 *
 * @param boozer Pointer to the data struct.
 * @param npsi Number of psi grid points in nu and theta data.
 * @param psi_min Minimum value in the psi grid.
 * @param psi_max Maximum value in the psi grid.
 * @param ntheta Number of boozer theta grid points in nu data.
 * @param nthetag Number of geometric theta grid points in theta data.
 * @param nu The difference between cylindrical angle phi and toroidal boozer.
 *           coordinate zeta, phi = zeta + nu [rad].
 * @param theta The boozer poloidal angle [rad].
 * @param nrzs The number of elements in rs and zs.
 * @param rs Separatrix contour R coordinates [m].
 * @param zs Separatrix contour z coordinates [m].
 *
 * @return Zero if initialization succeeded.
 */
int boozer_init(
    boozer_data *boozer, int npsi, real psi_min, real psi_max, int ntheta,
    int nthetag, int npadding, real *nu, real *theta, int nrzs, real *rs,
    real *zs);

/**
 * @brief Free allocated resources.
 *
 * Spline-interpolants are freed.
 *
 * @param boozer The struct whose fields are deallocated.
 */
void boozer_free(boozer_data *boozer);

/**
 * Offload data to the accelerator.
 *
 * @param boozer The struct to offload.
 */
void boozer_offload(boozer_data *boozer);

/**
 * Evaluate Boozer coordinates and partial derivatives.
 *
 * The output vector has the following elements:
 *
 * - psithetazeta[0]  = psi
 * - psithetazeta[1]  = dpsi/dR
 * - psithetazeta[2]  = dpsi/dphi
 * - psithetazeta[3]  = dpsi/dz
 * - psithetazeta[4]  = theta
 * - psithetazeta[5]  = dtheta/dR
 * - psithetazeta[6]  = dtheta/dphi
 * - psithetazeta[7]  = dtheta/dz
 * - psithetazeta[8]  = zeta
 * - psithetazeta[9]  = dzeta/dR
 * - psithetazeta[10] = dzeta/dphi
 * - psithetazeta[11] = dzeta/dz
 *
 * @param psithetazeta Evaluated Boozer coordinates and their gradients.
 * @param isinside A flag indicating whether the queried point was inside
 *        boozer grid.
 * @param r R coordinate.
 * @param phi phi coordinate.
 * @param z z coordinate.
 * @param bfield Pointer to magnetic field data.
 * @param boozer Pointer to boozerdata.
 *
 * @return Zero on success.
 */
DECLARE_TARGET_SIMD_UNIFORM(bfield, boozer)
a5err boozer_eval_psithetazeta(
    real psithetazeta[12], int *isinside, real r, real phi, real z,
    B_field_data *bfield, boozer_data *boozer);

#endif
