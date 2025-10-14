/**
 * @file boozer.h
 * Module for mapping cylindrical coordinates into straight-field-line (Boozer)
 * coordinates.
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "bfield.h"
#include "boozer_data.h"
#include "defines.h"
#include "utils/interp.h"

/**
 * Initialize the data for Boozer coordinate mapping.
 *
 * The tabulated values of the poloidal Boozer angle ``theta`` has a peculiarity
 * that ``theta`` is periodic but not in a typical sense: there is a
 * discontinuity when the geometrical angle crosses 2pi and ``theta`` drops from
 * 2pi to 0. The interpolant cannot deal with this, so we solve this by
 * requiring that the values of ``theta`` are padded in both directions. The
 * padded values do not include the "drop": instead they go above 2pi or below
 * 0. This way the splines are well behaved in the region [0, 2pi] and we never
 * interpolate outside this region. The padded values are counted in
 * ``nthetag``.
 *
 * @param boozer The struct to initialize.
 * @param npsi Number of psi grid points in ``nu`` and ``theta`` data.
 * @param ntheta Number of Boozer theta grid points in ``nu`` data.
 * @param nthetag Number of geometric theta grid points in ``theta`` data.
 * @param nrz The number of elements in ``rs`` and ``zs``.
 * @param npadding Number of padding values in ``theta`` in one direction.
 * @param psilim Limits of the uniform psi grid [Wb/rad].
 * @param nu The difference between cylindrical angle phi and toroidal boozer
 *        coordinate: nu = zeta - phi [rad].
 *        Layout: [psi_i, theta_j] = array[j*npsi + i].
 * @param theta The boozer poloidal angle [rad].
 *        Layout: [psi_i, thetag_j] = array[j*npsi + i].
 * @param rs Separatrix contour R coordinates [m].
 * @param zs Separatrix contour z coordinates [m].
 *
 * @return Zero if initialization succeeded.
 */
int Boozer_init(
    Boozer *boozer, size_t npsi, size_t ntheta, size_t nthetag, size_t nrz,
    int npadding, real psilim[2], real nu[npsi * ntheta],
    real theta[npsi * nthetag], real rs[nrz], real zs[nrz]);

/**
 * Free allocated resources.
 *
 * Spline-interpolants are freed.
 *
 * @param boozer The struct whose fields are deallocated.
 */
void Boozer_free(Boozer *boozer);

/**
 * Offload data to the accelerator.
 *
 * @param boozer The struct to offload.
 */
void Boozer_offload(Boozer *boozer);

DECLARE_TARGET_SIMD_UNIFORM(bfield, boozer)
/**
 * Find Boozer coordinates and their gradients at a given point in cylindrical
 * coordinates.
 *
 * @param psithetazeta Evaluated Boozer coordinates and their gradients
 *        [Wb, rad, rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz, theta, dtheta/dr,
 *        dtheta/dphi, dtheta/dz, zeta, dzeta/dr, dzeta/dphi, dzeta/dz].
 * @param isinside A flag indicating whether the queried point was inside
 *        the flux region where the Boozer coordinates are defined.
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param boozer The Boozer coordinate data.
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Boozer_map_coordinates(
    real psithetazeta[12], int *isinside, real r, real phi, real z, real t,
    Boozer *boozer, Bfield *bfield);

#endif
