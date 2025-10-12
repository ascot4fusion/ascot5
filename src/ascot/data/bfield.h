/**
 * @file bfield.h
 * Magnetic field interface.
 *
 * Provides functions and datatypes for magnetic field evaluation.
 */
#ifndef BFIELD_H
#define BFIELD_H

#include "bfield_data.h"
#include "defines.h"
#include "parallel.h"
#include "utils/interp.h"

/**
 * Free allocated resources.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void Bfield_free(Bfield *bfield);

/**
 * Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void Bfield_offload(Bfield *bfield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate poloidal flux.
 *
 * @param psi Evaluated poloidal flux [Wb/rad].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_psi(
    real psi[1], real r, real phi, real z, real t, Bfield *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate poloidal flux and its derivatives.
 *
 * @param psi_dpsi Evaluated poloidal flux and it's derivatives [Wb/rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, Bfield *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate normalized poloidal flux and its derivative (with respect to
 * poloidal flux) from poloidal flux.
 *
 * @param rho Evaluated normalized poloidal flux and it's derivatives [1].
 *        Layout: [rho, drho/dpsi].
 * @param psi Poloidal flux [Wb/rad].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_rho(real rho[2], real psi, Bfield *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate normalized poloidal flux and its derivatives from poloidal flux.
 *
 * @param rho_drho Evaluated normalized poloidal flux and it's derivatives [1].
 *        Layout: [rho, drho/dr, drho/dphi, drho/dz].
 * @param psi_dpsi Poloidal flux and it's derivatives [Wb/rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_rho_drho(
    real rho_drho[4], real psi_dpsi[4], Bfield *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate magnetic field vector.
 *
 * @param b Evaluated magnetic field vector [T].
 *        Layout: [br, bphi, bz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_b(
    real b[3], real r, real phi, real z, real t, Bfield *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate magnetic field vector and its derivatives.
 *
 * @param b_db Evaluated magnetic field vector and its derivatives [T].
 *        Layout: [br, dbr/dr, dbr/dphi, bz, dbz/dz, bphi, dbphi/dr, dbphi/dphi,
 *        dbphi/dz, bz, dbz/dr, dbz/dphi, dbz/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_b_db(
    real b_db[15], real r, real phi, real z, real t, Bfield *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate the magnetic axis (R, z) coordinates.
 *
 * Returns the position stored in the struct.
 *
 * @param axisrz Evaluated axis coordinates [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Bfield_eval_axis_rz(real rz[2], Bfield *bfield, real phi);
DECLARE_TARGET_END

#endif
