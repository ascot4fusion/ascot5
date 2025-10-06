/**
 * @file bfield_analytical.h
 * Analytical tokamak magnetic field implementation.
 *
 * This field combines the analytical equilibrium `[1]`_, which is quite good
 * model, to analytical ripple, which is significantly less realistic as it
 * is not divergence-free.
 *
 * .. _[1]: https://doi.org/10.1063/1.3328818
 */
#ifndef BFIELD_ANALYTICAL_H
#define BFIELD_ANALYTICAL_H
#include "defines.h"
#include "bfield.h"
#include "errors.h"
#include "parallel.h"

/**
 * Initialize the analytical magnetic field data.
 *
 * Assigns the fields in the struct with the provided values. All arrays in the
 * struct have fixed lengths so there's no need to allocate anything.
 *
 * @param bfield The struct to initialize.
 * @param nripple Number of tooridal field coils.
 * @param bphi Toroidal field strength at ``rmajor`` [T].
 * @param rmajor Plasma major radius [m].
 * @param rminor Plasma minor radius [m].
 * @param psiscaling  Scaling factor to scale poloidal flux [Wb/rad].
 * @param ripplescaling Ripple scaling factor.
 * @param rippledamping Ripple minor radius damping factor [m].
 * @param axisrz Magnetic axis (R, z) coordinates [m].
 * @param psilimits Poloidal flux at axis and separatrix [Wb/rad].
 * @param coefficients Coefficients for evaluating poloidal flux.
 *        Layout: [c1, c2, ..., c12, A].
 *
 * @return Zero if the initialization succeeded.
 */
int BfieldAnalytical_init(
    BfieldAnalytical *bfield, int nripple, real bphi, real rmajor, real rminor,
    real psiscaling, real ripplescaling, real rippledamping, real axisrz[2],
    real psilimits[2], real coefficients[13]);

/**
 * Free allocated resources.
 *
 * Does nothing since no resources were allocated during the initialization.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void BfieldAnalytical_free(BfieldAnalytical *bfield);

/**
 * Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void BfieldAnalytical_offload(BfieldAnalytical *bfield);

/**
 * Evaluate poloidal flux.
 *
 * @param psi Evaluated poloidal flux [Wb/rad].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldAnalytical_eval_psi(
    real psi[1], real r, real z, BfieldAnalytical *bfield);
DECLARE_TARGET_END

/**
 * Evaluate poloidal flux and its derivatives.
 *
 * @param psi_dpsi Evaluated poloidal flux and it's derivatives [Wb/rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldAnalytical_eval_psi_dpsi(
    real psi_dpsi[4], real r, real z, BfieldAnalytical *bfield);
DECLARE_TARGET_END

/**
 * Evaluate normalized poloidal flux and its derivatives.
 *
 * @param rho_drho Evaluated normalized poloidal flux and it's derivatives [1].
 *        Layout: [rho, drho/dr, drho/dphi, drho/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldAnalytical_eval_rho_drho(
    real rho_drho[4], real r, real z, BfieldAnalytical *bfield);
DECLARE_TARGET_END

/**
 * Evaluate magnetic field vector.
 *
 * If ``nripple`` is non-zero, the ripple contribution is included.
 *
 * @param b Evaluated magnetic field vector [T].
 *        Layout: [br, bphi, bz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldAnalytical_eval_b(
    real b[3], real r, real phi, real z, BfieldAnalytical *bfield);
DECLARE_TARGET_END

/**
 * Evaluate magnetic field vector and its derivatives.
 *
 * If ``nripple`` is non-zero, the ripple contribution is included.
 *
 * @param b_db Evaluated magnetic field vector and its derivatives [T].
 *        Layout: [br, dbr/dr, dbr/dphi, bz, dbz/dz, bphi, dbphi/dr, dbphi/dphi,
 *        dbphi/dz, bz, dbz/dr, dbz/dphi, dbz/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldAnalytical_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldAnalytical *bfield);
DECLARE_TARGET_END

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
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldAnalytical_eval_axisrz(real axisrz[2], BfieldAnalytical *bfield);
DECLARE_TARGET_END

#endif
