/**
 * @file bfield_spline2d.h
 * Spline-interpolated axisymmetric tokamak magnetic field implementation.
 *
 * Both poloidal flux and magnetic field is evaluated by interpolating
 * uniform 2D tables in which the data is provided. The poloidal component of
 * the evaluated magnetic field includes the contribution from the gradient of
 * the poloidal flux.
 */
#ifndef BFIELD_SPLINE2D_H
#define BFIELD_SPLINE2D_H
#include "bfield.h"
#include "defines.h"
#include "errors.h"
#include "parallel.h"

/**
 * Initialize the 2D spline-interpolated tokamak magnetic field data.
 *
 * Assigns the fields in the struct with the provided values, and initializes
 * and allocates the spline-interpolants.
 *
 * @param bfield The struct to initialize.
 * @param nr Number of R grid points in the data.
 * @param nz Number of z grid points in the data.
 * @param rlim Limits of the uniform R abscissa [m].
 * @param zlim Limits of the uniform z abscissa [m].
 * @param axisrz Magnetic axis (R, z) coordinates [m].
 * @param psilimits Poloidal flux at axis and separatrix [Wb/rad].
 * @param psi Tabulated values of poloidal flux [Wb/rad].
 *        Layout: (Ri, zj) = [j*nr + i] (C order).
 * @param br Tabulated values of R component of B [T].
 *        Layout: (Ri, zj) = [j*nr + i] (C order).
 * @param bz Tabulated values of z component of B [T].
 *        Layout: (Ri, zj) = [j*nr + i] (C order).
 * @param bphi Tabulated values of phi component of B [T].
 *        Layout: (Ri, zj) = [j*nr + i] (C order).
 *
 * @return  Zero if the initialization succeeded.
 */
int BfieldSpline2D_init(
    BfieldSpline2D *bfield, int nr, int nz, real rlim[2], real zlim[2],
    real axisrz[2], real psilimits[2], real psi[nr * nz], real br[nr * nz],
    real bz[nr * nz], real bphi[nr * nz]);

/**
 * Free allocated resources
 *
 * Spline-interpolants are freed.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void BfieldSpline2D_free(BfieldSpline2D *bfield);

/**
 * Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void BfieldSpline2D_offload(BfieldSpline2D *bfield);

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
a5err BfieldSpline2D_eval_psi(
    real psi[1], real r, real z, BfieldSpline2D *bfield);
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
a5err BfieldSpline2D_eval_psi_dpsi(
    real psi_dpsi[4], real r, real z, BfieldSpline2D *bfield);
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
a5err BfieldSpline2D_eval_rho_drho(
    real rho_drho[4], real r, real z, BfieldSpline2D *bfield);
DECLARE_TARGET_END

/**
 * Evaluate magnetic field vector.
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
a5err BfieldSpline2D_eval_b(
    real b[3], real r, real z, BfieldSpline2D *bfield);
DECLARE_TARGET_END

/**
 * Evaluate magnetic field vector and its derivatives.
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
a5err BfieldSpline2D_eval_b_db(
    real b_db[12], real r, real z, BfieldSpline2D *bfield);
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
a5err BfieldSpline2D_eval_axisrz(real axisrz[2], BfieldSpline2D *bfield);
DECLARE_TARGET_END
#endif
