/**
 * @file mhd_dynamic.h
 * Mhd perturbation where eigenmodes evolve in time.
 */
#ifndef MHD_DYNAMIC_H
#define MHD_DYNAMIC_H

#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd.h"
#include "parallel.h"

/**
 * Initialize dynamic MHD eigenmodes.
 *
 * @param mhd The struct to initialize.
 * @param nmode The number of modes.
 * @param nrho Number of points in rho grid.
 * @param moden Toroidal mode numbers.
 * @param modem Poloidal mode numbers.
 * @param rholim Minimum and maximum values in the uniform rho grid [1].
 * @param amplitude Mode amplitudes.
 * @param omega Mode frequencies [rad/s].
 * @param phase Mode phases [rad].
 * @param alpha Magnetic perturbation eigenfunctions [m].
 *        Layout: (rhoi, tj, mode) = [mode*nrho*ntime + j*nrho + i] (C order).
 * @param phi Electric perturbation eigenfunctions [V].
 *        Layout: (rhoi, tj, mode) = [mode*nrho*ntime + j*nrho + i] (C order).
 *
 * @return Zero if the initialization succeeded.
 */
int MhdDynamic_init(
    MhdDynamic *mhd, size_t n, size_t nrho, size_t ntime, int nmode[n],
    int mmode[n], real rholim[2], real tlim[2], real amplitude[n],
    real omega[n], real phase[n], real alpha[n * nrho * ntime],
    real phi[n * nrho * ntime]);

/**
 * Free allocated resources.
 *
 * @param mhd The struct whose fields are deallocated.
 */
void MhdDynamic_free(MhdDynamic *mhd);

/**
 * Offload data to the accelerator.
 *
 * @param mhd The struct to offload.
 */
void MhdDynamic_offload(MhdDynamic *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd, bfield, boozer, include_mode)
/**
 * Evaluate the MHD terms used in the guiding center equations of motion.
 *
 * @param alpha The magnetic perturbation term and it's derivatives [m].
 * @param Phi The electric perturbation term and it's derivatives [V].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param include_mode Include perturbation from just this mode number.
 *        Default: MHD_INCLUDE_ALL.
 * @param mhd The MHD data.
 * @param bfield The magnetic field data.
 * @param boozer The Boozer coordinate data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t MhdDynamic_eval_alpha_Phi(
    real alpha[5], real Phi[5], real r, real phi, real z, real t,
    size_t include_mode, MhdDynamic *mhd, Bfield *bfield, Boozer *boozer);

#endif
