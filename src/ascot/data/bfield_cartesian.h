/**
 * @file bfield_cartesian.h
 * Cartesian magnetic field implementation.
 *
 * This is a magnetic field that is defined on a Cartesian basis by a vector
 * at the origo and a Jacobian. This field doesn't have flux surfaces (as this
 * is mainly for testing) so the axis location and the poloidal flux value
 * (which is constant in space) are arbitrary.
 */
#ifndef BFIELD_CARTESIAN_H
#define BFIELD_CARTESIAN_H
#include "bfield.h"
#include "defines.h"
#include "parallel.h"

/**
 * Initialize the Cartesian magnetic field.
 *
 * Assign the fields in the struct with the provided values. All arrays in the
 * struct have fixed lengths so there's no need to allocate anything.
 *
 * @param bfield The struct to initialize.
 * @param psival Value that is returned when quering the magnetic flux [Vs/m].
 * @param rhoval Value that is returned when quering the normalized poloidal
 *        flux [1].
 * @param axisrz (R, z) coordinates that are returned when the magnetic axis
 *        location is queried [m].
 * @param bxyz Magnetic field vector at origo [T].
 *        Layout: [bx, by, bz].
 * @param jacobian Magnetic field Jacobian [T/m]
 *        Layout: [dbx/dx, dbx/dy, dbx/dz, dby/dx, dby/dy, dby/dz, dbz/dx,
 *        dbz/dy, dbz/dz].
 *
 * @return Zero if the initialization succeeded.
 */
int BfieldCartesian_init(
    BfieldCartesian *bfield, real psival, real rhoval, real axisrz[2],
    real bxyz[3], real jacobian[9]);

/**
 * Free allocated resources.
 *
 * Does nothing since no resources were allocated during the initialization.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void BfieldCartesian_free(BfieldCartesian *bfield);

/**
 * Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void BfieldCartesian_offload(BfieldCartesian *bfield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate poloidal flux.
 *
 * The poloidal flux is a constant in space.
 *
 * @param psi Evaluated poloidal flux [Wb/rad].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldCartesian_eval_psi(real psi[1], BfieldCartesian *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate poloidal flux and its derivatives.
 *
 * The poloidal flux is a constant in space and it's derivatives are zero.
 *
 * @param psi_dpsi Evaluated poloidal flux and it's derivatives [Wb/rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldCartesian_eval_psi_dpsi(real psi_dpsi[4], BfieldCartesian *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate magnetic field vector.
 *
 * @param b Evaluated magnetic field vector in cylindrical basis [T].
 *        Layout: [br, bphi, bz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldCartesian_eval_b(
    real b[3], real r, real phi, real z, BfieldCartesian *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate magnetic field vector and its derivatives.
 *
 * @param b_db Evaluated magnetic field vector and its derivatives in
 *        cylindrical basis [T].
 *        Layout: [br, dbr/dr, dbr/dphi, bz, dbz/dz, bphi, dbphi/dr, dbphi/dphi,
 *        dbphi/dz, bz, dbz/dr, dbz/dphi, dbz/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldCartesian_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldCartesian *bfield);
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
err_t BfieldCartesian_eval_axisrz(real axisrz[2], BfieldCartesian *bfield);
DECLARE_TARGET_END

#endif
