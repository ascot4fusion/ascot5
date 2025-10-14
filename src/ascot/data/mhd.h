/**
 * @file mhd.h
 * MHD perturbation interface.
 *
 * Provides functions and data types for evaluating MHD eigenmodes.
 */
#ifndef MHD_H
#define MHD_H

#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd_data.h"

/**
 * ``includemode`` parameter to include all modes (default)
 */
#define MHD_INCLUDE_ALL 0

/**
 * Free allocated resources.
 *
 * @param mhd The struct whose fields are deallocated.
 */
void Mhd_free(Mhd *mhd);

/**
 * Offload data to the accelerator.
 *
 * @param mhd The struct to offload.
 */
void Mhd_offload(Mhd *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd, bfield, boozer, include_mode)
/**
 * Evaluate the MHD terms used in the guiding center equations of motion.
 *
 * @param alpha The magnetic perturbation term and it's derivatives [m].
 *        Layout: [alpha, dalpha/dt, dalpha/dr, dalpha/dphi, dalpha/dz].
 * @param Phi The electric perturbation term and it's derivatives [V].
 *        Layout: [Phi, dPhi/dt, dPhi/dr, dPhi/dphi, dPhi/dz].
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
err_t Mhd_eval_alpha_Phi(
    real alpha[5], real Phi[5], real r, real phi, real z, real t,
    size_t include_mode, Mhd *mhd, Bfield *bfield, Boozer *boozer);

DECLARE_TARGET_SIMD_UNIFORM(mhd, bfield, boozer, pertonly, include_mode)
/**
 * Evaluate perturbed fields and electric potential explicitly.
 *
 * @param b The magnetic field perturbation [T].
 *        Layout: [br, bphi, bz].
 * @param e The electric field perturbation [V/m].
 *        Layout: [er, ephi, ez].
 * @param Phi The electric potential perturbation [V].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param include_background If true, include the background field in ``b``.
 * @param include_mode Include perturbation from just this mode number.
 *        Default: MHD_INCLUDE_ALL.
 * @param mhd Mhd data.
 * @param bfield Magnetic field data.
 * @param boozer The Boozer coordinate data.
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
err_t Mhd_eval_perturbation(
    real b[3], real e[3], real Phi[1], real r, real phi, real z, real t,
    int include_background, size_t include_mode, Mhd *mhd, Bfield *bfield,
    Boozer *boozer);

DECLARE_TARGET_SIMD_UNIFORM(mhd)
/**
 * Get number of modes.
 *
 * @param mhd MHD data.
 *
 * @return Number of modes.
 */
size_t Mhd_get_n_modes(Mhd *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd)
/**
 * Get mode toroidal numbers.
 *
 * @param mhd MHD data.
 *
 * @return Mode toroidal numbers.
 */
const int *Mhd_get_nmode(Mhd *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd)
/**
 * Get mode poloidal numbers.
 *
 * @param mhd MHD data.
 *
 * @return Mode poloidal numbers.
 */
const int *Mhd_get_mmode(Mhd *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd)
/**
 * Get mode amplitudes.
 *
 * @param mhd MHD data.
 *
 * @return Mode amplitudes [1].
 */
const real *Mhd_get_amplitude(Mhd *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd)
/**
 * Get mode frequencies.
 *
 * @param mhd MHD data.
 *
 * @return Mode frequencies [rad/s].
 */
const real *Mhd_get_frequency(Mhd *mhd);

DECLARE_TARGET_SIMD_UNIFORM(mhd)
/**
 * Get mode phases.
 *
 * @param mhd MHD data.
 *
 * @return Mode phases [rad].
 */
const real *Mhd_get_phase(Mhd *mhd);
#endif
